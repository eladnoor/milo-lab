#!/usr/bin/python

"""\
html_writer.py - Construct HTML pages

"""

import os
import sys

class HtmlWriter:
    def __init__(self, filename, force_path_creation=True, flush_always=True):
        self.filename = filename
        self.filepath = os.path.dirname(filename)
        self.flush_always = flush_always
        if (not os.path.exists(self.filepath)):
            if (force_path_creation and not os.path.exists(self.filepath)):
                os.mkdir(self.filepath)
            else:
                raise Exception("cannot write to HTML file %s since the directory doesn't exist" % filename)
        
        self.file = open(self.filename, "w")
        self.write('<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">\n')
        self.write('<head>\n')
        self.write('<script type="text/javascript" src="expandCollapse.js"></script>\n')
        self.write('</head>\n')
        self.write('<html>\n<body>\n')
        self.write_js()
    
    def write_js(self):
        if (os.path.exists(self.filepath + '/expandCollapse.js')):
            return
        file = open(self.filepath + '/expandCollapse.js', 'w')
        file.write("""function toggleMe(a){
  var e=document.getElementById(a);
  if(!e)return true;
  if(e.style.display=="none"){
    e.style.display="block"
  } else {
    e.style.display="none"
  }
  return true;
}
""")
        file.close()
    
    def write(self, str):
        if (self.file == None):
            raise Exception("cannot write to this HTML since it is already closed")
        self.file.write(str)
        if (self.flush_always):
            self.flush()
    
    def write_ol(self, l):
        self.write("<ol>\n")
        for mem in l:
            self.write("  <li>%s</li>\n" % mem)
        self.write("</ol>\n")
            
    def write_ul(self, l):
        self.write("<ul>\n")
        for mem in l:
            self.write("  <li>%s</li>\n" % mem)
        self.write("</ul>\n")
    
    def flush(self):
        if (self.file != None):
            self.file.flush()
    
    def embed_img(self, fig_fname, alternative_string=""):
        self.write('<img src="' + fig_fname + '" atl="' + alternative_string + '" />')
    
    def embed_svg(self, fig_fname, alternative_string="", width=320, height=240, name=''):
        self.write('<object data="%s" type="image/svg+xml" width="%d" height="%d" name="%s"/></object>' % (fig_fname, width, height, name))
    
    def branch(self, relative_path, link_text=None):
        """
            Branches the HTML file by creating a new HTML and adding a link to it with the desired text
        """
        if (link_text == None):
            link_text = relative_path
            
        self.write("<a href=\"" + relative_path + ".html\">" + link_text + "</a>")
        return HtmlWriter(os.path.join(self.filepath, relative_path + ".html"))
    
    def close(self):
        self.write("</body>\n</html>\n")
        self.file.flush()
        self.file.close()
        self.file = None

def test():
    html_write = HtmlWriter("../res/test.html")
    html_write.write("hello world")

if __name__ == '__main__': test()