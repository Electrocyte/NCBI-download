import os
import time
import json
import socket
import argparse
from ftplib import FTP, error_perm, error_temp

# https://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/Xylariaceae_sp._FL1272/all_assembly_versions/GCA_022432385.1_XyFL1272_2/

def establish_ftp(timeout=120):  # 120 seconds timeout
    ftp = FTP('ftp.ncbi.nlm.nih.gov', timeout=timeout)
    ftp.login()
    ftp.cwd(base_url)
    return ftp


def retry_on_errors(func):
    def wrapper(*args, **kwargs):
        MAX_ATTEMPTS = 5  # for example, adjust as needed
        RETRY_DELAY = 5  # delay between attempts in seconds, adjust as needed
        attempts = 0
        while True:
            try:
                return func(*args, **kwargs)
            except (EOFError, BrokenPipeError, socket.timeout, error_temp) as e:  # add socket.timeout here
                attempts += 1
                if attempts > MAX_ATTEMPTS:
                    print(f"Maximum retry attempts reached. Failed to execute {func.__name__}.")
                    raise  # re-raise the last exception

                print(f"Caught exception: {str(e)}, attempting to reestablish connection (Attempt {attempts}).")  # print the exception
                global ftp  # make sure to use the global ftp instance

                ftp = establish_ftp()  # reestablish the connection
                ftp.cwd(base_url)  # navigate back to the base URL

                print(f"Waiting {RETRY_DELAY} seconds before next attempt...")
                time.sleep(RETRY_DELAY)  # delay before the next attempt
    return wrapper


@retry_on_errors
def cwd(ftp_instance, *args):
    ftp_instance.cwd(*args)

@retry_on_errors
def nlst(ftp_instance, *args):
    return ftp_instance.nlst(*args)

@retry_on_errors
def retrbinary(ftp_instance, *args, **kwargs):
    ftp_instance.retrbinary(*args, **kwargs)



parser = argparse.ArgumentParser(description='Download gbff files from NCBI')

parser.add_argument('-t', '--target', help='Input kingdom target from NCBI')
parser.add_argument('-o', '--output', help='directory top level')
parser.add_argument('-f', '--failures', help='JSON file with previous failures')
parser.add_argument('-s', '--skip', action='store_true', help='Skip existing folders')

args = parser.parse_args()

target = args.target
output = args.output
failures_file = args.failures
skip = args.skip

output_dir = f"{output}/{target}"
base_url = f"/genomes/genbank/{target}"

ftp = establish_ftp()
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

if failures_file is not None:
    with open(failures_file, 'r') as f:
        failures_dict = json.load(f)

    # Delete only failed files listed in failures file
    for failure in failures_dict.keys():
        if os.path.exists(failure):
            os.remove(failure)

    directories = []
    for failure in failures_dict.keys():
        failure = failure.replace(output, "")
        failure = failure.replace(target, "")
        failure = failure.split('/')[1]
        if len(failure) > 0:
            directories.append(failure)

else:
    directories = nlst(ftp)  # use the decorated function

print(f"First 10 items: {directories[0:10]}; Total: {len(directories)}")

MAX_RETRIES = 4

n = 0
failures = 0
failure_list = []
seen_genus_species = set()
for directory in directories:

    if skip:
        if os.path.exists(os.path.join(output_dir, directory)):
            n += 1
            print(f"Skipping {directory}")
            continue

    # Extract genus species
    genus_species = '_'.join(directory.split('_')[:2]) # Change this according to the exact format of your directories

    # If we've seen this genus species before, skip this iteration
    if genus_species in seen_genus_species:
        continue
    seen_genus_species.add(genus_species)

    # Save current directory
    current_directory = ftp.pwd()
    print(f"Processing {directory} ({n+1} of {len(directories)})")

    directory_to_change = directory + '/all_assembly_versions'
    print(f"Changing to {directory_to_change}")
    try:
        cwd(ftp, directory_to_change)
    except error_perm:
        print(f"Could not change to {directory_to_change}, skipping.")
        continue
    time.sleep(0.05)  # add a delay between FTP requests

    subdirectories = nlst(ftp)  # use the decorated function

    # for subdirectory in subdirectories:
    for subdirectory in subdirectories[:1]:  # Only process the first subdirectory
        while True:  # Loop to attempt reconnection if FTP server disconnects
            try:
                current_subdirectory = ftp.pwd()
                break  # If this line is reached, no error was raised, so exit the loop
            except (EOFError, error_temp) as e:
                print(f"Caught exception: {e}, attempting to reestablish connection.")
                ftp = establish_ftp()
            ftp.cwd(current_directory + '/' + directory_to_change + '/' + subdirectory)

        # Check if local files exist
        taxon_dir = os.path.join(output_dir, directory)
        version_dir = os.path.join(taxon_dir, subdirectory)
        if failures_file is not None:
            if os.path.exists(version_dir):
                # Assuming a directory is considered complete if it exists
                print(f"Local files in {version_dir} already exists, skipping.")
                continue

        try:
            ftp.cwd(subdirectory)
        except error_perm:
            print(f"Could not change to {subdirectory}, skipping.")
            # Change back to previous directory
            ftp.cwd(current_subdirectory)
            continue

        try:
            filenames = nlst(ftp)  # use the decorated function
        except (EOFError, BrokenPipeError, socket.timeout) as e:
            print(f"Maximum retry attempts reached for listing files in directory {directory}. Moving on to next directory.")
            continue

        for filename in filenames:
            if filename.endswith('genomic.gbff.gz'):

                taxon_dir = os.path.join(output_dir, directory)
                if not os.path.exists(taxon_dir):
                    os.makedirs(taxon_dir)

                version_dir = os.path.join(taxon_dir, subdirectory)
                if not os.path.exists(version_dir):
                    os.makedirs(version_dir)

                local_filename = os.path.join(version_dir, filename)
                if os.path.exists(local_filename or local_filename.replace(".gz", "")):  # Skip downloading if the file already exists
                    print(f"{filename} already exists, skipping.")
                    continue

                print(f"Downloading {filename} to {local_filename}")
                retries = 0
                while retries < MAX_RETRIES:
                    try:
                        with open(local_filename, 'wb') as file_handle:
                            retrbinary(ftp, 'RETR ' + filename, file_handle.write)
                        break  # success, move on to next file
                    except Exception as e:
                        retries += 1
                        print(f"Caught exception: {e}, attempting to reestablish connection (Attempt {retries}).")
                        failures += 1
                        failure_list.append(filename)
                        ftp = establish_ftp()
                        ftp.cwd(base_url + '/' + directory + '/all_assembly_versions/' + subdirectory)
                        if retries >= MAX_RETRIES:
                            print(f"Failed to download {filename} after {MAX_RETRIES} attempts, skipping.")
                            break

                print(f"Output file size: {os.stat(local_filename).st_size/1e6:.1e} MB")

        # When done with subdirectory, change back to previous directory
        ftp.cwd(current_subdirectory)
    n += 1
    # Change back to previous directory
    ftp.cwd(current_directory)
ftp.quit()

print(f"Downloaded {n} organisms with {failures} failures.")
print(seen_genus_species)

save_fails = f'{output_dir}/failures.txt'
with open(save_fails, 'w') as file_handle:
    print(f"Saving failures to {save_fails}")
    file_handle.write('\n'.join(failure_list))

# time python ./download-genbank.py
