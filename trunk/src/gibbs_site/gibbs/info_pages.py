from django.shortcuts import render_to_response


def AboutPage(request):
    """Renders the landing page."""
    return render_to_response('about.html', {})


def FAQPage(request):
    """Renders the landing page."""
    return render_to_response('faq.html', {})