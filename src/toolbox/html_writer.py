#!/usr/bin/python

"""\
html_writer.py - Construct HTML pages

"""

import os
import xml.dom.minidom

class BaseHtmlWriter:
    def __init__(self):
        pass

    def write(self):
        raise Exception("class not implemented")

    def write_header(self):
        self.write('<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">\n')
        self.write('<head>\n')
        self.write('<script type="text/javascript" src="expandCollapse.js"></script>\n')
        self.write('</head>\n')
        self.write('<html>\n<body>\n')

    def write_js(self, path):
        if (os.path.exists(path + '/expandCollapse.js')):
            return
        file = open(path + '/expandCollapse.js', 'w')
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
        
    def write_ol(self, l):
        self.write("<ol>\n")
        for mem in l:
            self.write("  <li>%s</li>\n" % str(mem))
        self.write("</ol>\n")
            
    def write_ul(self, l):
        self.write("<ul>\n")
        for mem in l:
            self.write("  <li>%s</li>\n" % str(mem))
        self.write("</ul>\n")
    
    def embed_img(self, fig_fname, alternative_string=""):
        self.write('<img src="' + fig_fname + '" atl="' + alternative_string + '" />')
    
    def embed_svg(self, fig_fname, alternative_string="", width=320, height=240, name=''):
        self.write('<object data="%s" type="image/svg+xml" width="%d" height="%d" name="%s"/></object>' % (fig_fname, width, height, name))
    
    def embed_matplotlib_figure(self, fig, width=320, height=240):
        """
            Converts the figure to an SVG DOM and uses the inline SVG option to 
            add it directly into the HTML (without creating a separate SVG file).
        """
        fig.savefig('.svg', format='svg')
        self.extract_svg_from_file('.svg', width, height)
        os.remove('.svg')

    def embed_dot(self, Gdot, width=320, height=240):
        """
            Converts the DOT graph to an SVG DOM and uses the inline SVG option to 
            add it directly into the HTML (without creating a separate SVG file).
        """
        Gdot.write('.svg', prog='dot', format='svg')
        self.extract_svg_from_file('.svg', width, height)
        os.remove('.svg')
        
    def extract_svg_from_file(self, fname, width=320, height=240):
        x = xml.dom.minidom.parse(fname)
        svg = x.getElementsByTagName("svg")[0]
        svg.setAttribute('width', '%dpt' % width)
        svg.setAttribute('height', '%dpt' % height)
        self.write(svg.toprettyxml(indent='  ', newl=''))
    
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

class NullHtmlWriter(BaseHtmlWriter):
    def __init__(self):
        BaseHtmlWriter.__init__(self)
    
    def write(self, str):
        pass

class HtmlWriter(BaseHtmlWriter):
    def __init__(self, filename, force_path_creation=True, flush_always=True):
        BaseHtmlWriter.__init__(self)
        self.filename = filename
        self.filepath = os.path.dirname(filename)
        self.flush_always = flush_always
        if (not os.path.exists(self.filepath)):
            if (force_path_creation and not os.path.exists(self.filepath)):
                os.mkdir(self.filepath)
            else:
                raise Exception("cannot write to HTML file %s since the directory doesn't exist" % filename)
        
        self.file = open(self.filename, "w")
        self.write_header()
        self.write_js(self.filepath)
    
    def write(self, str):
        if (self.file == None):
            raise Exception("cannot write to this HTML since it is already closed")
        self.file.write(str)
        if (self.flush_always):
            self.file.flush()
            
    def close(self):
        BaseHtmlWriter.close(self)
        self.file.flush()
        self.file.close()
        self.file = None

def test():
    html_write = HtmlWriter("../res/test.html")
    html_write.write("<h1>hello world</h1>\n")
    import pylab
    fig = pylab.figure()
    pylab.plot([1, 2, 3, 4], [4, 3, 2, 1], 'g--')
    html_write.embed_matplotlib_figure(fig, 1000, 1000)

if __name__ == '__main__': test()