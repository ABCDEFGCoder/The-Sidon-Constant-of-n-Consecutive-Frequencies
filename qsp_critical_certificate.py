"""Compatibility interface for the printed degree-seven uniqueness box.

The global bridge is rechecked from the unrestricted phase-cover source.
The extra narrow RUR evaluation box is proved to lie in the printed box.
"""
from flint import fmpq,ctx
from critical_bridge import check_printed_box

def verify(c,bridge=True):
    assert c['format']=='Sn_uniform_QSP_critical_root_v1' and c['n']==7
    box,bounds=check_printed_box(7,global_bridge=bridge)
    radius=fmpq(c['bridge_radius'])
    for pair,X in zip(c['centre'],box):
        centre=fmpq(*pair)
        assert X.lower()<centre-radius and centre+radius<X.upper()
    return {'global_bridge_checked':bool(bridge),'printed_box_bounds':bounds}
