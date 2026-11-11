"""Additive ISV circRNA evidence overlay; never assigns biological validation."""
from math import isfinite

RED = '#c62828'


def marker_specs(molecule, lo, hi):
    """Return donor/acceptor markers. Coordinates are 1-based inclusive bases.

    No endpoint snapping, inferential filtering, or evidence-status promotion.
    A cropped endpoint is omitted from the viewport only, not the evidence.
    """
    specs = []
    for joint in molecule.get('circular_junctions', []):
        if joint.get('coordinate_system') != '1-based-inclusive':
            raise ValueError('Circular junction requires explicit 1-based-inclusive coordinates')
        status = joint.get('validation_status', 'candidate')
        if status not in ('candidate', 'orthogonally_validated'):
            raise ValueError('Unrecognized circular validation status')
        if status == 'orthogonally_validated' and not joint.get('validation_source'):
            raise ValueError('Orthogonal validation requires a source record')
        strand = joint.get('strand')
        if strand not in ('+', '-'):
            raise ValueError('Circular junction requires strand')
        donor, acceptor = joint['donor'], joint['acceptor']
        for position in (donor, acceptor):
            if isinstance(position, bool) or not isinstance(position, int) or position < 1:
                raise ValueError('Circular endpoint must be a positive integer')
        if not (donor > acceptor if strand == '+' else donor < acceptor):
            raise ValueError('Circular endpoint order disagrees with transcript strand')
        for role, position in [('donor', donor), ('acceptor', acceptor)]:
            if lo <= position <= hi:
                specs.append(dict(role=role, position=position, validation_status=status,
                                  linestyle='solid'))
    return specs


def draw_circular_markers(ax, molecule, xz, y_center, exon_height, lo, hi):
    """Hollow red circles whose outer diameter equals exon height in display units."""
    from matplotlib.patches import Circle
    if not isfinite(exon_height) or exon_height <= 0:
        raise ValueError('Invalid exon height')
    artists = []
    for spec in marker_specs(molecule, lo, hi):
        artist = Circle((xz(spec['position']), y_center), radius=exon_height / 2,
                        facecolor='white', edgecolor=RED,
                        linewidth=min(1.2, exon_height / 5), linestyle=spec['linestyle'],
                        zorder=8, clip_on=True)
        artist.set_gid('circular_' + spec['role'] + '_' + str(spec['position']))
        ax.add_patch(artist)
        artists.append(artist)
    return artists


def validate_circular_payload(res):
    """Avoid the native renderer's silent row cap for additive circular views."""
    for panel in res.get('panels', []):
        molecules = panel.get('molecules', [])
        if any(m.get('circular_junctions') for m in molecules):
            if len(molecules) > 4000:
                raise ValueError('Circular view exceeds 4000 rows; paginate explicitly, preserving evidence')
            if res.get('has_multiple_loci'):
                raise ValueError('Circular overlay expects a single genomic coordinate frame')
            for molecule in molecules:
                marker_specs(molecule, 1, float('inf'))
