from django.conf.urls.defaults import *
from django.conf import settings

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('',
    (r'^media/(?P<path>.*)$', 'django.views.static.serve', {'document_root': settings.MEDIA_ROOT }),
    (r'^$', 'gibbs_site.gibbs.main_page.MainPage'),
    (r'^about', 'gibbs_site.gibbs.info_pages.AboutPage'),
    (r'^faq', 'gibbs_site.gibbs.info_pages.FAQPage'),
    (r'^classic_reactions', 'gibbs_site.gibbs.info_pages.ClassicReactions'),
    (r'^compound', 'gibbs_site.gibbs.compound_page.CompoundPage'),
    (r'^enzyme', 'gibbs_site.gibbs.enzyme_page.EnzymePage'),
    (r'^reaction', 'gibbs_site.gibbs.reaction_page.ReactionPage'),
    (r'^graph_reaction', 'gibbs_site.gibbs.reaction_graph.ReactionGraph'),
    (r'^search', 'gibbs_site.gibbs.search_results_page.ResultsPage'),
    (r'^suggest', 'gibbs_site.gibbs.suggest.SuggestJson'),
    (r'^robots\.txt', 'gibbs_site.gibbs.info_pages.Robots'),
    # Example:
    # (r'^gibbs_site/', include('gibbs_site.foo.urls')),

    # Uncomment the admin/doc line below and add 'django.contrib.admindocs' 
    # to INSTALLED_APPS to enable admin documentation:
    # (r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    (r'^admin/', include(admin.site.urls)),
)
