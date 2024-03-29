import os
from datetime import datetime

import click
import yaml
import logging

from vrs_anvil import Manifest, run_command_in_background, get_process_info
from vrs_anvil.annotator import annotate_all
from logging.handlers import RotatingFileHandler
import pathlib

# Set up logging
log_format = "%(asctime)s %(threadName)s %(name)s [%(levelname)s] %(message)s"

_logger = logging.getLogger("vrs_anvil.cli")


@click.group(invoke_without_command=True)
@click.version_option(package_name="vrs_anvil")
@click.option(
    "--manifest",
    type=click.Path(exists=False),
    default="manifest.yaml",
    help="Path to manifest file.",
)
@click.option("--verbose", default=False, help="Log more information", is_flag=True)
@click.option("--max_errors", default=10, help="Number of acceptable errors.")
@click.pass_context
def cli(ctx, verbose: bool, manifest: str, max_errors: int):
    """GA4GH GKS utility for AnVIL."""

    _log_level = logging.INFO
    if verbose:
        _log_level = logging.DEBUG

    timestamp_str = datetime.now().strftime("%Y%m%d_%H%M%S")

    try:
        with open(manifest, "r") as stream:
            manifest = Manifest.model_validate(yaml.safe_load(stream))

            # only create a persistent log for annotate subcommand
            if ctx.invoked_subcommand == "annotate":

                # Create a rotating file handler with a max size of 10MB and keep 3 backup files
                log_path = pathlib.Path(manifest.state_directory) / f"vrs_anvil_{timestamp_str}_{os.getpid()}.log"
                file_handler = RotatingFileHandler(
                    log_path, maxBytes=10 * 1024 * 1024, backupCount=3
                )
                file_handler.setLevel(_log_level)
                file_handler.setFormatter(logging.Formatter(log_format))

                # Add the file handler to the logger
                # logger.addHandler(file_handler)
                # basicConfig call removed, which prevents the default configuration that logs to the console.
                logging.basicConfig(
                    level=_log_level, format=log_format, handlers=[file_handler]
                )

                click.secho(
                    f"🪵  Logging to {log_path}, level {logging.getLevelName(_log_level)}",
                    fg="yellow",
                )

            else:
                logging.basicConfig(level=_log_level, format=log_format)

            ctx.ensure_object(dict)
            ctx.obj["manifest"] = manifest
            ctx.obj["verbose"] = verbose
            ctx.obj["max_errors"] = max_errors
            ctx.obj["timestamp_str"] = timestamp_str

            if verbose:
                click.secho(f"📢  {manifest}", fg="green")

    except Exception as exc:
        click.secho(f"{exc}", fg="yellow")
        ctx.ensure_object(dict)


@cli.command("annotate")
@click.option('--scatter', help='Start a background process per VCF file.', required=False, default=False,
              is_flag=True, show_default=True)
@click.pass_context
def annotate_cli(ctx, scatter: bool):
    """Read manifest file, annotate variants, all parameters controlled by manifest.yaml."""

    timestamp_str = ctx.obj["timestamp_str"]

    if not scatter:
        try:
            assert "manifest" in ctx.obj, "Manifest not found."
            manifest = ctx.obj["manifest"]
            _logger.debug(f"Manifest: {ctx.obj['manifest']}")
            click.secho("🚧  annotating variants", fg="yellow")
            metrics_file = annotate_all(manifest, max_errors=ctx.obj["max_errors"], timestamp_str=timestamp_str)
            click.secho(f"📊  metrics available in {metrics_file}", fg="green")
        except Exception as exc:
            click.secho(f"{exc}", fg="red")
            _logger.exception(exc)
    else:
        try:
            assert "manifest" in ctx.obj, "Manifest not found."
            parent_manifest = ctx.obj["manifest"]
            c = 0
            scattered_processes = []
            for _ in parent_manifest.vcf_files:
                # create a new manifest for each VCF file, based on the parent manifest, clone the parent manifest
                child_manifest = Manifest.parse_obj(parent_manifest.model_dump())
                child_manifest.vcf_files = [_]
                child_manifest.num_threads = 1
                child_manifest.disable_progress_bars = True
                child_manifest_path = pathlib.Path(child_manifest.work_directory) / f"manifest_{timestamp_str}_{c}.yaml"
                with open(child_manifest_path, "w") as stream:
                    yaml.dump(child_manifest.model_dump(), stream)
                child_pid = run_command_in_background(f'vrs_anvil --manifest {child_manifest_path} annotate')
                click.secho(f"🚧  annotating {_} on pid {child_pid}", fg="yellow")
                scattered_processes.append({'pid': child_pid, 'manifest': str(child_manifest_path), 'vcf': _})
                c += 1
            scattered_processes_path = pathlib.Path(parent_manifest.work_directory) / f"scattered_processes_{timestamp_str}.yaml"
            scattered_processes = {'parent_pid': os.getpid(), 'processes': scattered_processes}
            with open(scattered_processes_path, "w") as stream:
                yaml.dump(scattered_processes, stream)
            click.secho(f"📊 scattered processes available in {scattered_processes_path}", fg="green")
        except Exception as exc:
            click.secho(f"{exc}", fg="red")
            _logger.exception(exc)


@cli.command("ps")
@click.pass_context
def ps_cli(ctx):
    """Show status of latest scatter command."""
    try:
        assert "manifest" in ctx.obj, "Manifest not found."
        parent_manifest = ctx.obj["manifest"]
        scattered_processes_path = pathlib.Path(parent_manifest.work_directory)
        scattered_processes_paths = sorted(x for x in scattered_processes_path.glob("scattered_processes_*.yaml"))
        if not scattered_processes_paths:
            click.secho(f"🚧  no scattered processes found in {parent_manifest.work_directory}/scattered_processes_*.yaml", fg="red")
            return
        scattered_processes_path = scattered_processes_paths[-1]
        state_dir = pathlib.Path(parent_manifest.state_directory)
        with open(scattered_processes_path, "r") as stream:
            scattered_processes = yaml.safe_load(stream)
            for _ in scattered_processes['processes']:
                log_file = 'NA'
                metrics_file = 'NA'
                try:
                    log_file = sorted(state_dir.glob(f"vrs_anvil*{_['pid']}.log"))[-1]
                    metrics_file = sorted(state_dir.glob(f"metrics_*{_['pid']}.yaml"))[-1]
                except IndexError:
                    pass

                click.secho(f"🚧  pid: {str(_['pid'])}, manifest: {str(_['manifest'])}, vcf: {str(_['vcf'])}, metrics_file: {metrics_file}, log_file: {log_file}", fg="yellow")
                process = get_process_info(_['pid'])
                if not process:
                    # assert pathlib.Path(metrics_file).exists(), f"metrics file not found: {metrics_file}"
                    # assert pathlib.Path(log_file).exists(), f"log file not found: {log_file}"
                    click.secho("  ✅  completed", fg="green")
                else:
                    io_counters = 'NA'
                    memory_info = 'NA'
                    if process.status() == 'running':
                        try:
                            if hasattr(process, 'io_counters'):
                                io_counters = process.io_counters()
                            if hasattr(process, 'memory_info'):
                                memory_info = process.memory_info()
                        except Exception as exc:
                            _logger.info(f"could not get io_counters/memory_info pid: {_['pid']} error:{exc}")
                    click.secho(f"  📊 {process.status()} cpu_percent: {process.cpu_percent(interval=0.1)}%, io_counters: {io_counters}, memory_info: {memory_info}", fg="yellow")
    except Exception as exc:
        click.secho(f"{exc}", fg="red")
        _logger.exception(exc)


if __name__ == "__main__":
    cli()
