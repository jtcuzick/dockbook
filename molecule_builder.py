

import streamlit.components.v1 as components

def draw_molecule():
    components.html("""
        <!DOCTYPE html>
        <html>
        <head>
            <script src="https://unpkg.com/kekule/dist/kekule.js?modules=chemWidget"></script>
            <link href="https://unpkg.com/kekule/dist/themes/default/kekule.css" rel="stylesheet">
            <style>
                #sketcher {
                    width: 100%;
                    height: 400px;
                    border: 1px solid #ddd;
                    margin-bottom: 10px;
                }
            </style>
        </head>
        <body>
            <div id="sketcher"></div>
            <script>
                var editor = new Kekule.Editor.StructureEditor(document.getElementById('sketcher'));
                editor.setChemObj(new Kekule.Molecule());

                window.getSmiles = function() {
                    var mol = editor.getChemObj();
                    return Kekule.IO.saveFormatData(mol, 'smi');
                };
            </script>
        </body>
        </html>
    """, height=420)

