import logging
import pylab
import PIL, PIL.Image, StringIO

from django.http import Http404
from django.http import HttpResponse
from django.shortcuts import render_to_response
from gibbs import concentration_profile
from gibbs import reaction
from gibbs import reaction_form

from matplotlib import pyplot
import pylab
from mpl_toolkits.mplot3d import Axes3D
import random
import itertools


def ReactionGraph(request):    
    """Renders a page for a particular reaction."""
    form = reaction_form.ReactionForm(request.GET)
    if not form.is_valid():
        logging.error(form.errors)
        raise Http404

    rxn = reaction.Reaction.FromForm(form)
    
    fig = pylab.figure()
    ax = Axes3D(fig)
    phs = pylab.arange(0.001, 10.0, 0.25)
    pmgs = pylab.arange(0.001, 10.0, 0.25)
    pairs = itertools.product(phs, pmgs)
    pts = [(p, m, rxn.DeltaGTag(pH=p, pMg=m)) for p,m in pairs]
    phs = [a[0] for a in pts]
    pmgs = [a[1] for a in pts]
    dgs = [a[2] for a in pts]
    ax.scatter(phs, pmgs, dgs)
    ax.set_xlabel('pH')
    ax.set_ylabel('pMg')
    ax.set_zlabel('dG\'')
    
    buffer = StringIO.StringIO()
    canvas = pyplot.get_current_fig_manager().canvas
    canvas.draw()
    pilImage = PIL.Image.fromstring("RGB", canvas.get_width_height(), canvas.tostring_rgb())
    pilImage.save(buffer, "PNG")
    pylab.close()
    
    return HttpResponse(buffer.getvalue(), mimetype="image/png")
