import click

seed = click.option('-S', '--seed', type=int, default=None,
              help='Seed for random number generator (e.g. 17, or 2017-04-11).'
              ' Defaults to system time.')
              
feature_dim = click.option('-d','--feature_dim', type=click.IntRange(0, 21),
    default=1,
    help='Max dependent interaction order/dimension considered when constructing position features.')

directions = click.option('--directions',
    type=click.Choice(['AtoC', 'AtoG', 'AtoT', 'CtoA', 'CtoG', 'CtoT',
                       'GtoA', 'GtoC', 'GtoT', 'TtoA', 'TtoC', 'TtoG',
                       'All', 'None']), default='All',
    help='Examples of choices of mut_directions option are:\
     - All for all directions, \
     - specific direction (e.g. AtoC) for a specified direction, or \
     - None for no direction type spedified when building a classifier')
