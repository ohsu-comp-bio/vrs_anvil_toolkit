
import click
import yaml
import logging

from vrs_anvil import Manifest

_logger = logging.getLogger(__name__)


@click.group(invoke_without_command=True)
@click.version_option(package_name='vrs_anvil')
@click.option('--manifest_path', type=click.Path(exists=True), default='manifest.yaml', help='Path to manifest file.')
@click.option('--verbose', default=False, help='Log more information', is_flag=True)
@click.pass_context
def cli(ctx, verbose: bool, manifest_path: str):
    """GA4GH GKS utility for AnVIL."""
    if verbose:
        logging.basicConfig(level=logging.DEBUG)
        _logger.setLevel(logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
        _logger.setLevel(logging.INFO)

    try:
        with open(manifest_path, 'r') as stream:
            manifest = Manifest.parse_obj(yaml.safe_load(stream))
    except Exception as exc:
        click.secho(f"Error in initializing: {exc}", fg='red')
        _logger.exception(exc)
        exit(1)

    if not ctx.invoked_subcommand:
        click.secho(manifest, fg='green')
        click.secho('No subcommand invoked', fg='yellow')

    ctx.ensure_object(dict)
    ctx.obj['manifest'] = manifest
    ctx.obj['verbose'] = verbose


@cli.command('annotate')
@click.argument('metrics_path', default='metrics.tsv')
@click.pass_context
def annotate_cli(ctx, metrics_path: str):
    """Read manifest file, annotate variants

    \b
    metrics_path: path to output file.
    """
    try:
        _logger.debug(f"Manifest: {ctx.obj['manifest']}")
        click.secho(f"TODO - 🚧 annotate variants, write results to {metrics_path}", fg='yellow')
    except Exception as exc:
        click.secho(f"Error in processing: {exc}", fg='red')
        _logger.exception(exc)


if __name__ == '__main__':
    cli()
