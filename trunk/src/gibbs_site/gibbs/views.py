import logging
from django.shortcuts import render_to_response
from django.http import Http404
from gibbs import compound_form
from gibbs import concentration_profile
from gibbs import models
from gibbs import reaction_form
from gibbs import search_form
from gibbs import service_config

def GetBalancednessWarning(reaction, ph=None, ionic_strength=None,
                           concentration_profile=None):
    if reaction.IsBalanced():
        return None
    
    warning = 'Reaction is not balanced!'
    if not reaction.CanBalanceWithWater():
        return warning
    
    url = reaction.GetBalanceWithWaterLink(ph, ionic_strength,
                                           concentration_profile)
    warning += ' <a href="%s">Balance with water?</url>' % url
    return warning


def CompoundPage(request):
    """Renders a page for a particular compound."""
    form = compound_form.CompoundForm(request.GET)
    if not form.is_valid():
        logging.error(form.errors)
        raise Http404
    
    # Compute the delta G estimate.
    kegg_id = form.cleaned_compoundId
    compound = models.Compound.objects.get(kegg_id=kegg_id)
    compound.StashTransformedSpeciesEnergies(form.cleaned_ph,
                                             form.cleaned_ionic_strength)
    delta_g_estimate = compound.DeltaG(
        pH=form.cleaned_ph, ionic_strength=form.cleaned_ionic_strength)
    
    template_data = {'compound': compound, 
                     'ph': form.cleaned_ph,
                     'ionic_strength': form.cleaned_ionic_strength,
                     'delta_g_estimate': delta_g_estimate,
                     'no_dg_explanation': compound.no_dg_explanation,
                     'kegg_link': compound.GetKeggLink()}
    return render_to_response('compound_page.html', template_data)


def ReactionPage(request):
    """Renders a page for a particular reaction."""
    form = reaction_form.ReactionForm(request.GET)
    if not form.is_valid():
        logging.error(form.errors)
        raise Http404

    i_s = form.cleaned_ionic_strength
    ph = form.cleaned_ph
    
    clean_reactants = form.cleaned_reactantIds
    clean_products = form.cleaned_productIds
    
    cprofile_name = form.cleaned_concentration_profile
    cprofile = concentration_profile.GetProfile(cprofile_name)
    
    reactant_names = form.cleaned_reactantNames
    product_names = form.cleaned_productNames
    
    zipped_reactants = zip(form.cleaned_reactantCoeffs, clean_reactants, reactant_names)
    zipped_products = zip(form.cleaned_productCoeffs, clean_products, product_names)
    
    reaction = models.Reaction.FromIds(zipped_reactants, zipped_products)
    if form.cleaned_balance_w_water:
        reaction.TryBalanceWithWater()
    
    # Compute the delta G estimate.
    delta_g_estimate = reaction.DeltaG(pH=ph,
                                       ionic_strength=i_s,
                                       concentration_profile=cprofile)
    
    warning = GetBalancednessWarning(reaction, ph, i_s, cprofile_name)
    template_data = {'reaction': reaction,
                     'query': form.cleaned_query,
                     'ph': form.cleaned_ph,
                     'ionic_strength': form.cleaned_ionic_strength,
                     'delta_g_estimate': delta_g_estimate,
                     'concentration_profile': cprofile_name,
                     'warning': warning}
    return render_to_response('reaction_page.html', template_data)

    
def ResultsPage(request):
    """Renders the search results page for a given query."""
    form = search_form.SearchForm(request.GET)
    if not form.is_valid():
        raise Http404
    
    query_parser = service_config.Get().query_parser
    reaction_matcher = service_config.Get().reaction_matcher
    matcher = service_config.Get().compound_matcher
    
    query = str(form.cleaned_query)
    ph = form.cleaned_ph
    ionic_strength = form.cleaned_ionic_strength
    template_data = {'query': query,
                     'ph': ph,
                     'ionic_strength': ionic_strength}
    
    # Check if we should parse and process the input as a reaction.
    if query_parser.IsReactionQuery(query):
        try:
            parsed_reaction = query_parser.ParseReactionQuery(query)
        except Exception:
            template_data['warning'] = 'Failed to parse your search query!'
            return render_to_response('error_page.html', template_data)

        reaction_matches = reaction_matcher.MatchReaction(parsed_reaction)
        best_reaction = reaction_matches.GetBestMatch()
        
        if not best_reaction:
            logging.error('Failed to match reaction query.')
            raise Http404
        
        reactants, products = best_reaction
        
        reaction = models.Reaction.FromIds(reactants, products)
        delta_g_estimate = reaction.DeltaG(
            pH=ph, ionic_strength=ionic_strength)

        template_data['warning'] = GetBalancednessWarning(reaction)  
        template_data.update({'delta_g_estimate': delta_g_estimate,
                              'reaction': reaction})
        return render_to_response('reaction_page.html', template_data)

    else:
        # Otherwise we try to parse it as a single compound.
        results = matcher.Match(query)
        delta_g_estimate = None
        compound = results[0].value
        compound.StashTransformedSpeciesEnergies(ph, ionic_strength)
        if results:
            delta_g_estimate = compound.DeltaG(
                pH=ph, ionic_strength=ionic_strength)
        
        template_data['kegg_link'] = results[0].value.GetKeggLink()
        template_data['results'] = results
        template_data['delta_g_estimate'] = delta_g_estimate
        return render_to_response('search_results.html',
                                  template_data)

    raise Http404
