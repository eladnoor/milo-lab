#!/usr/bin/python

from util import django_utils

# NOTE(flamholz): This is crappy. We're using the real database for
# a unit test. I wish I knew of a better way.
django_utils.SetupDjango()

from gibbs import models

def main():
    urls = ['http://milolab.webfactional.com/\n',
            'http://milolab.webfactional.com/faq\n',
            'http://milolab.webfactional.com/about\n']
    
    for compound in models.Compound.objects.all():
        urls.append('http://milolab.webfactional.com%s\n' % compound.link)
        
    for reaction in models.StoredReaction.objects.all():
        urls.append('http://milolab.webfactional.com%s\n' % reaction.link)
    
    for enzyme in models.Enzyme.objects.all():
        urls.append('http://milolab.webfactional.com%s\n' % enzyme.link)
    
    
    f = open('sitemap.txt', 'w')
    f.writelines(urls)
    f.close()
            

if __name__ == '__main__':
    main()
                
                
