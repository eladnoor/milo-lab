from django.conf.urls.defaults import patterns, url
import primers.views
import settings

# Uncomment the next two lines to enable the admin:
# from django.contrib import admin
# admin.autodiscover()

urlpatterns = patterns('',
    # Examples:
    # url(r'^$', 'mysite.views.home', name='home'),
    # url(r'^mysite/', include('mysite.foo.urls')),

    # Uncomment the admin/doc line below to enable admin documentation:
    # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    # url(r'^admin/', include(admin.site.urls)),
    url(r'^media/(?P<path>.*)$', 'django.views.static.serve', {
            'document_root': settings.MEDIA_ROOT, }),
    ('^$', primers.views.welcome),
    ('^primero/search/$', primers.views.primer_search),
    ('^primero/query/$', primers.views.results),
    ('^primero/update_loc/$', primers.views.update_location),
    ('^primero/fetch_login/$', primers.views.fetch_login),
    ('^primero/fetchmail/$', primers.views.fetch_mail),
    ('^primero/pending/$', primers.views.pending),
    ('^primero/update/$', primers.views.update_pending),
    ('^primero/error/$', primers.views.example_error),
    ('^primero/upload/$', primers.views.upload_sigma_pdf),
)
