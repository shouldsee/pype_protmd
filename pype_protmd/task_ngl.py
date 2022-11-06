# import os
# __file__ = os.path.realpath(__file__)

from .depend_mol import know_ngl,Controller
ctl = Controller()
know_ngl(ctl)
ctl.build()

from jinja2 import Template
import jinja2

#loadFile = "rcsb://1crn"
# loadFile = 'https://files.rcsb.org/download/1PGB.pdb'
# loadFile = '/pype/build/1PGB.pdb'
loadFile = '/pype/build/1PGB_solv_ions.gro'
def minimal_ngl_html(loadFile):
    HEADER_JS = ''
    HEADER_JS += f'''
    loadFile = {loadFile!r}
    '''

    HEADER_JS += '''

    var stage;
    document.addEventListener("DOMContentLoaded", function (){
        var stage = new NGL.Stage( "viewport" );

        window.addEventListener( "resize", function( event ){
            stage.handleResize();
        }, false );

        stage.loadFile( loadFile, { defaultRepresentation: true } );
    })
    '''

    HEADER_CSS = '''
    * { margin: 0; padding: 0; }
    html, body { width: 100%; height: 100%; overflow: hidden; }
    '''

    #'<script src="./node_modules/ngl/dist/ngl.js"></script>'
    #PATH_NGL = './node_modules/ngl/dist'
    buf = f'''
    <html>
    <head>
    <script src="{ctl['init_ngl'].check_ctx}"></script>
    <script>
    {HEADER_JS}
    </script>

    <style>
    {HEADER_CSS}
    </style>
    </head>
    <body>
    <div id="viewport" style="width:100%; height:100%;"></div>
    </body># os.getcwd())
    # assert 0, NGL_REL

    </html>
    '''
import os

with open(__file__ + '.template.jinja2') as file_:
    TEMPLATE_NGL  = Template(file_.read(), undefined=jinja2.StrictUndefined)

def get_ngl_html(loadFile, TARGET_FILE = None):
    if TARGET_FILE is None:
        TARGET_FILE  =os.path.basename(__file__)+'.html'
    TARGET_DIR = os.path.dirname(os.path.realpath(TARGET_FILE))

    NGL_REL = os.path.relpath( ctl.state['git/ngl'].built+'/examples', TARGET_DIR)
    buf = TEMPLATE_NGL.render(NGL_REL=NGL_REL,loadFile=loadFile)
    with open(TARGET_FILE, 'w') as f:
        f.write(buf)
    return TARGET_FILE


if __name__=='__main__':
    TARGET_FILE = get_ngl_html(loadFile)
    TARGET_FILE_REL = os.path.relpath(TARGET_FILE,'.')
    with open('index.html','w') as f:
        f.write(f'''
    <a href="{TARGET_FILE_REL}">{TARGET_FILE_REL}</a>
        ''')
    # print(buf)
