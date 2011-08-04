from django.shortcuts import render_to_response


def AboutPage(request):
    """Renders the about page."""
    return render_to_response('about.html', {})


def FAQPage(request):
    """Renders the FAQ page."""
    return render_to_response('faq.html', {})


def RedoxReview(request):
    """Renders robots.txt."""
    return render_to_response('redox_review.html', {})


def ClassicReactions(request):
    """Renders the classic reactions page."""
    return render_to_response('classic_reactions.html', {})


def DownloadPage(request):
    """Renders the download page."""
    
    ph_values = [5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0]
    ph_values = map(lambda x: '%.1f' % x, ph_values)
    return render_to_response('download.html', {'ph_values': ph_values})
    

def Robots(request):
    """Renders robots.txt."""
    return render_to_response('robots.txt', {})