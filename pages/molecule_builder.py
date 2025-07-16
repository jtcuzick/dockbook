import streamlit as st
import streamlit.components.v1 as components

def draw_molecule():
    components.html("""
        <html>
        <head>
            <script src="https://unpkg.com/kekule@1.0.2/dist/kekule.min.js"></script>
            <link rel="stylesheet" href="https://unpkg.com/kekule@1.0.2/dist/themes/default/kekule.css" />
        </head>
        <body>
            <div id="chemDraw" style="width: 100%; height: 500px;"></div>
            <script>
                var editor = new Kekule.Editor.Composer(document.getElementById('chemDraw'));
                editor.setPredefinedSetting('fullFunc');
            </script>
        </body>
        </html>
    """, height=520)


