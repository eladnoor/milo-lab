#!/usr/bin/python

import os.path

from django import template
from django.conf import settings
from toolbox import util

settings.configure()

def render_to_string(template_name, data=None):
    """Renders a template of the given name to a string.
    
    Args:
        template_name: the name of a template file in pygibbs/templates.
        data: a dictionary of template data.
    
    Returns:
        The string of the rendered template.
    """
    template_filename = os.path.abspath('pygibbs/templates/' + template_name)
    my_data = data or {}
    
    c = template.Context(my_data)
    t = template.Template(open(template_filename).read())
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