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
    flow = commands.add_parser("catalog-flow", help="Linear/chimeric catalog flow diagram from a box/value table")
    from .catalog_flow import add_arguments as add_flow_arguments
    add_flow_arguments(flow)
    flow.set_defaults(func=run_catalog_flow)
    ideo = commands.add_parser("chromosome-ideogram",
                               help="Banded GRCh38 ideograms (SVG) of a trans-splicing junction and its best DNA contact")
    from .chromosome_ideogram import add_arguments as add_ideogram_arguments
    add_ideogram_arguments(ideo)
    ideo.set_defaults(func=run_ideogram)


def run_isv(args):
    from .isv_web.isv_web_remote import run_from_args
    return run_from_args(args)


def run_psi(args):
    from .psi_boxplot import run
    return run(args)


def run_psi_grid(args):
    from .psi_boxplot_grid import run
    return run(args)


def run_ideogram(args):
    from .chromosome_ideogram import run
    return run(args)


def run_catalog_flow(args):
    from .catalog_flow import run
    return run(args)
