from django.shortcuts import render_to_response


def AboutPage(request):
    """Renders the about page."""
    return render_to_response('about.html', {})


def FAQPage(request):
    """Renders the FAQ page."""
    return render_to_response('faq.html', {})


def ClassicReactions(request):
    """Renders the classic reactions page."""
    return render_to_response('classic_reactions.html', {})


def DownloadPage(request):
    """Renders the download page."""
    return render_to_response('download.html', {})
    

def Robots(request):
    """Renders robots.txt."""
    return render_to_response('robots.txt', {})