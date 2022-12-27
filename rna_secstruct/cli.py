import click 


"""
# multi commmand format
@click.group()
def cli():
    pass


@cli.command()
def func1():
    pass
"""

@click.command()
def cli():
    pass


if __name__ == '__main__':
    cli()
