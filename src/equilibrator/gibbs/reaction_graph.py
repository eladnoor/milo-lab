import json
import logging

from gibbs import reaction
from gibbs import reaction_graph_form
from django.shortcuts import render_to_response


def ReactionGraph(request):
    """Renders the graph page."""
    form = reaction_graph_form.ReactionGraphForm(request.GET)
    if not form.is_valid():
        logging.error(form.errors)
        raise Http404

    mode = 'varyPh'
    if form.cleaned_vary_is:
        mode = 'varyIs'

    rxn = reaction.Reaction.FromForm(form)
    reactant_data = json.dumps([c.ToJson() for c in rxn.filtered_reactants])
    product_data = json.dumps([c.ToJson() for c in rxn.filtered_products])
    template_data = {"reactant_data": reactant_data,
                     "product_data": product_data,
                     "mode": mode,
                     "reaction": rxn}
    
    return render_to_response('reaction_graph.html', template_data)