import click


@click.group()
@click.version_option(package_name='vrs_anvil')
@click.pass_context
def cli(ctx):
    pass


@cli.command()
def to_do():
    """To do."""
    pass


if __name__ == '__main__':
    cli()
