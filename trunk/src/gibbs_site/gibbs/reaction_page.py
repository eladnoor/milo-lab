import logging

from django.http import Http404
from django.shortcuts import render_to_response
from gibbs import concentration_profile
from gibbs import reaction
from gibbs import reaction_form


_REACTION_TEMPLATES_BY_SUBMIT = {'Update': 'reaction_page.html',
                                 'Save': 'print_reaction.html'}


def ReactionPage(request):    
    """Renders a page for a particular reaction."""
    form = reaction_form.ReactionForm(request.GET)
    if not form.is_valid():
        logging.error(form.errors)
        raise Http404
    
    # Figure out which template to render (based on which submit button they pressed).
    template_name = _REACTION_TEMPLATES_BY_SUBMIT.get(form.cleaned_submit,
                                                      'reaction_page.html')

    # Fetch parameters and compounds.
    i_s = form.cleaned_ionic_strength
    ph = form.cleaned_ph
    
    clean_reactants = form.cleaned_reactantIds
    clean_products = form.cleaned_productIds
    all_ids = clean_reactants + clean_products
    
    # Fetch custom concentrations if any.
    reactant_concentrations = form.cleaned_reactantConcentrations
    product_concentrations = form.cleaned_productConcentrations
    all_concentrations = reactant_concentrations + product_concentrations
    
    # Build the appropriate concentration profile.
    cprofile_name = form.cleaned_concentration_profile
    cprofile = concentration_profile.GetProfile(
        cprofile_name, all_ids, all_concentrations)
    
    reactant_names = form.cleaned_reactantNames
    product_names = form.cleaned_productNames
    
    # Build the Reaction object and potentially balance it.
    zipped_reactants = zip(form.cleaned_reactantCoeffs, clean_reactants, reactant_names)
    zipped_products = zip(form.cleaned_productCoeffs, clean_products, product_names)
    rxn = reaction.Reaction.FromIds(zipped_reactants, zipped_products,
                                    concentration_profile=cprofile,
                                    pH=ph, ionic_strength=i_s)
    query = form.cleaned_query
    if form.cleaned_balance_w_water:
        rxn.TryBalanceWithWater()
        query = rxn.GetQueryString()
    if form.cleaned_balance_electrons:
        rxn.BalanceElectrons()
        query = rxn.GetQueryString()
    
    # Render the template.
    balance_with_water_link = rxn.GetBalanceWithWaterLink(query)
    balance_electrons_link = rxn.GetBalanceElectronsLink(query)
    template_data = {'reaction': rxn,
                     'query': query,
                     'ph': form.cleaned_ph,
                     'ionic_strength': form.cleaned_ionic_strength,
                     'concentration_profile': cprofile_name,
                     'balance_with_water_link': balance_with_water_link,
                     'balance_electrons_link': balance_electrons_link}
    return render_to_response(template_name, template_data)