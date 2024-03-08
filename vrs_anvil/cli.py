
import click
import yaml
import logging

from vrs_anvil import Manifest
from vrs_anvil.annotator import annotate_all
from logging.handlers import RotatingFileHandler
import pathlib

# Set up logging
log_format = "%(asctime)s %(filename)s [%(levelname)s] %(message)s"


_logger = logging.getLogger(__name__)


@click.group(invoke_without_command=True)
@click.version_option(package_name='vrs_anvil')
@click.option('--manifest', type=click.Path(exists=False), default='manifest.yaml', help='Path to manifest file.')
@click.option('--verbose', default=False, help='Log more information', is_flag=True)
@click.option('--max_errors', default=10, help='Number of acceptable errors.')
@click.pass_context
def cli(ctx, verbose: bool, manifest: str, max_errors: int):
    """GA4GH GKS utility for AnVIL."""

    _log_level = logging.INFO
    if verbose:
        _log_level = logging.DEBUG

    try:
        with open(manifest, 'r') as stream:
            manifest = Manifest.parse_obj(yaml.safe_load(stream))
            # basicConfig call removed, which prevents the default configuration that logs to the console.
            # logging.basicConfig(level=_log_level, format=log_format)

            # Create a rotating file handler with a max size of 10MB and keep 3 backup files
            log_path = pathlib.Path(manifest.state_directory) / "vrs_anvil.log"
            file_handler = RotatingFileHandler(
                    log_path, maxBytes=10 * 1024 * 1024, backupCount=3)
            file_handler.setLevel(_log_level)
            file_handler.setFormatter(logging.Formatter(log_format))

            # Add the file handler to the logger
            logger = logging.getLogger()
            logger.addHandler(file_handler)
            click.secho(f"🪵  Logging to {log_path}, level {logging.getLevelName(_log_level)}", fg='yellow')

            ctx.ensure_object(dict)
            ctx.obj['manifest'] = manifest
            ctx.obj['verbose'] = verbose
            ctx.obj['max_errors'] = max_errors

            if verbose:
                click.secho(f"📢  {manifest}", fg='green')

    except Exception as exc:
        click.secho(f"{exc}", fg='yellow')
        ctx.ensure_object(dict)


@cli.command('annotate')
@click.pass_context
def annotate_cli(ctx):
    """Read manifest file, annotate variants, all parameters controlled by manifest.yaml."""

    try:
        assert 'manifest' in ctx.obj, "Manifest not found."
        manifest = ctx.obj['manifest']
        _logger.debug(f"Manifest: {ctx.obj['manifest']}")
        click.secho("🚧  annotating variants", fg='yellow')
        metrics_file = annotate_all(manifest, max_errors=ctx.obj['max_errors'])
        click.secho(f"🥳  metrics available in {metrics_file}", fg='green')
    except Exception as exc:
        click.secho(f"{exc}", fg='red')
        _logger.exception(exc)


if __name__ == '__main__':
    cli()
