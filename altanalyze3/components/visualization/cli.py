"""AltAnalyze3 visualization commands; no alternative read-selection algorithms."""


def add_arguments(parser):
    from .isv_web.isv_web_remote import build_arg_parser
    from .psi_boxplot import add_arguments as add_psi_arguments

    commands = parser.add_subparsers(dest="visualization_action", required=True)
    isv = commands.add_parser("isv-domains", parents=[build_arg_parser()], add_help=False,
                              help="Existing native ISV renderer, including saved-input replay")
    isv.set_defaults(func=run_isv)
    psi = commands.add_parser("psi-boxplot", help="Donor PSI boxes from complete native isoform count tables")
    add_psi_arguments(psi)
    psi.set_defaults(func=run_psi)
    grid = commands.add_parser("psi-boxplot-grid", help="One-page panels of literal native donor PSI")
    from .psi_boxplot_grid import add_arguments as add_grid_arguments
    add_grid_arguments(grid)
    grid.set_defaults(func=run_psi_grid)


def run_isv(args):
    from .isv_web.isv_web_remote import run_from_args
    return run_from_args(args)


def run_psi(args):
    from .psi_boxplot import run
    return run(args)


def run_psi_grid(args):
    from .psi_boxplot_grid import run
    return run(args)
