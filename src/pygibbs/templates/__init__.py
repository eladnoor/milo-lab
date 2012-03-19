#!/usr/bin/python

from django.conf import settings

settings.configure(TEMPLATE_DIRS=('pygibbs/templates',))

import os.path

from django.template import loader, Context
from toolbox import util



def render_to_string(template_name, data=None):
    """Renders a template of the given name to a string.
    
    Args:
        template_name: the name of a template file in pygibbs/templates.
        data: a dictionary of template data.
    
    Returns:
        The string of the rendered template.
    """
    my_data = data or {}
    c = Context(my_data)
    t = loader.get_template(template_name)
    return t.render(c)

def render_to_file(template_name, data, output_filename):
    """Renders a template to a given file.
    
    Will create the parent directory of the output file if not present.
    
    Args:
        template_name: the name of a template file in pygibbs/templates.
        data: a dictionary of template data.
        output_filename: the name/path of the file to write to.
    """
    dir = os.path.abspath(os.path.dirname(output_filename))
    if not os.path.exists(dir):
        util._mkdir(dir)
        
    open(output_filename, 'w').write(render_to_string(
        template_name, data))