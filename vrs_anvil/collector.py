import os
from typing import Generator
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed
from vrs_anvil import Manifest
from google.cloud import storage
import boto3


def download_s3_object(bucket_name, object_name, destination_file_name) -> str:
    """Download an object from an S3 bucket and save it to a file."""
    if os.path.exists(destination_file_name):
        return destination_file_name
    s3 = boto3.client('s3')
    s3.download_http_file(bucket_name, object_name, destination_file_name)
    return destination_file_name


def download_google_blob(bucket_name, source_blob_name, destination_file_name) -> str:
    """Downloads a blob from the bucket."""
    if os.path.exists(destination_file_name):
        return destination_file_name

    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)
    blob = bucket.blob(source_blob_name)
    blob.download_to_filename(destination_file_name)
    return destination_file_name


def download_http_file(url, destination_dir) -> str:
    """Download a file from a URL and save it to a directory."""
    filename = os.path.join(destination_dir, url.split("/")[-1])
    # check if file already exists
    if os.path.exists(filename):
        return filename
    response = requests.get(url)
    with open(filename, 'wb') as f:
        f.write(response.content)
    return filename


def create_symlink_to_work_directory(work_directory: str, vcf_file: str) -> str:
    """Create a link to the work directory for vcf file."""
    # Extract the filename from the vcf_file path

    vcf_filename = os.path.basename(vcf_file)

    # Construct the full path for the symlink
    symlink_path = os.path.join(work_directory, vcf_filename)

    # Do not re-create it if it already exists
    if os.path.islink(symlink_path):
        return symlink_path

    # Create the symlink
    vcf_file = os.path.abspath(vcf_file)
    os.symlink(vcf_file, symlink_path)

    return symlink_path


def collect_manifest_urls(manifest: Manifest) -> Generator[str, None, None]:
    """Collect the URLs from the manifest and download them."""
    with ThreadPoolExecutor(max_workers=manifest.num_threads) as executor:
        futures = []
        for vcf_file in manifest.vcf_files:
            if vcf_file.startswith("http"):
                futures.append(executor.submit(download_http_file, vcf_file, manifest.work_directory))
            elif vcf_file.startswith("gs://"):
                bucket_name, blob_name = vcf_file[5:].split("/", 1)
                destination_file_name = os.path.join(manifest.work_directory, blob_name)
                futures.append(executor.submit(download_google_blob, bucket_name, blob_name, destination_file_name))
            elif vcf_file.startswith("s3://"):
                bucket_name, object_name = vcf_file[5:].split("/", 1)
                destination_file_name = os.path.join(manifest.work_directory, object_name)
                futures.append(executor.submit(download_s3_object, bucket_name, object_name, destination_file_name))
            elif vcf_file.startswith("file://"):
                futures.append(executor.submit(create_symlink_to_work_directory, manifest.work_directory, vcf_file[7:]))
            else:
                futures.append(executor.submit(create_symlink_to_work_directory, manifest.work_directory, vcf_file))
        # Yield results as they become available
        for _ in as_completed(futures):
            yield _.result()
