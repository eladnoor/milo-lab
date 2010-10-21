from django.conf.urls.defaults import *
from django.conf import settings

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('',
    (r'^media/(?P<path>.*)$', 'django.views.static.serve', {'document_root': settings.MEDIA_ROOT }),
    (r'^$', 'gibbs_site.gibbs.main_page.MainPage'),
    (r'^compound', 'gibbs_site.gibbs.views.CompoundPage'),
    (r'^reaction', 'gibbs_site.gibbs.views.ReactionPage'),
    (r'^search', 'gibbs_site.gibbs.views.ResultsPage'),
    (r'^suggest', 'gibbs_site.gibbs.suggest.SuggestJson'),
    # Example:
    # (r'^gibbs_site/', include('gibbs_site.foo.urls')),

    # Uncomment the admin/doc line below and add 'django.contrib.admindocs' 
    # to INSTALLED_APPS to enable admin documentation:
    # (r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    (r'^admin/', include(admin.site.urls)),
)
