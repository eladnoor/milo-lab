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
        rxn = reaction.Reaction.FromIds(reactants, products,
                                        concentration_profile=cprofile,
                                        pH=ph, ionic_strength=ionic_strength)
        
        balance_with_water_link = rxn.GetBalanceWithWaterLink(query)
        balance_electrons_link = rxn.GetBalanceElectronsLink(query)
        template_data.update({'reaction': rxn,
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
