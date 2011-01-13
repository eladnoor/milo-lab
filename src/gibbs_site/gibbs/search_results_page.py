import logging

from django.http import Http404
from django.shortcuts import render_to_response
from gibbs import concentration_profile
from gibbs import reaction
from gibbs import search_form
from gibbs import service_config

    
def ResultsPage(request):
    """Renders the search results page for a given query."""
    form = search_form.SearchForm(request.GET)
    if not form.is_valid():
        raise Http404
    
    query_parser = service_config.Get().query_parser
    reaction_matcher = service_config.Get().reaction_matcher
    matcher = service_config.Get().compound_matcher
    
    query = form.cleaned_query
    logging.info(query)
    logging.info('Is reaction query: %s' % query_parser.IsReactionQuery(query))
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
            return render_to_response('error_page.html', template_data)

        reaction_matches = reaction_matcher.MatchReaction(parsed_reaction)
        best_reaction = reaction_matches.GetBestMatch()
        
        if not best_reaction:
            logging.error('Failed to match reaction query.')
            raise Http404
        
        reactants, products = best_reaction
        
        cprofile = concentration_profile.GetProfile()
        rxn = reaction.Reaction.FromIds(reactants, products, cprofile)
        delta_g_estimate = rxn.DeltaG(
            pH=ph, ionic_strength=ionic_strength)

        balance_with_water_link = rxn.GetBalanceWithWaterLink(
            ph, ionic_strength, cprofile.name, query)
        balance_electrons_link = rxn.GetBalanceElectronsLink(
            ph, ionic_strength, cprofile.name, query)
        template_data.update({'delta_g_estimate': delta_g_estimate,
                              'no_dg_explanation': rxn.NoDeltaGExplanation(),
                              'reaction': rxn,
                              'balance_with_water_link': balance_with_water_link,
                              'balance_electrons_link': balance_electrons_link})
        return render_to_response('reaction_page.html', template_data)

    else:
        # Otherwise we try to parse it as a single compound.
        results = matcher.Match(query)        
        template_data['results'] = results
        return render_to_response('search_results.html',
                                  template_data)

    raise Http404
