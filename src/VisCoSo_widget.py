import anywidget
import traitlets
import ipywidgets as widgets
from IPython.display import display

from .cola_parser import Parser
from .problem import EmptyException
from .util import *


def viscoso():
    Widget()


class Widget(object):
    def __init__(self):

        self.program = None
        self.text_area = None

        text = widgets.HTML(
            """
            <script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>
            <script type="text/x-mathjax-config">
                MathJax.Hub.Config({TeX: {extensions: ["action.js"] }});
            </script>
            """
        )
        display(text)

        self.input_layout()
        self.output_layout()

    def input_layout(self):

        self.text_area = widgets.Textarea(
            value="A CoLa program",
            placeholder="Type a CoLa program",
            description="",
            disabled=False,
        )
        msg = widgets.Label("Or write a CoLa program in the textbox:")
        ub = self.upload_button()
        input_layout = widgets.VBox([ub, msg, self.text_area])
        display(input_layout)

    def output_layout(self):
        solve_button = widgets.Button(
            description="Solve",
            disabled=False,
            button_style="",
            tooltip="Solve",
            icon="check",
        )

        display(solve_button)
        output = widgets.Output()

        def on_button_solve(b):
            if self.program is None:
                self.program = self.text_area.value

            parser = Parser(cola=self.program)
            try:
                parser.parse()
                sol = parser.problem.solve(debug=False)
                self.viscoso_widget(sol)
            except EmptyException:
                with output:
                    print("Could not find a problem :(")
                display(output)

        solve_button.on_click(on_button_solve)

    def upload_button(self):
        upload_button = widgets.FileUpload(
            accept="",  # Accepted file extension e.g. '.txt', '.pdf', 'image/*', 'image/*,.pdf'
            multiple=False,  # True to accept multiple files upload else False
        )
        msg = widgets.Label("Upload a CoLa program:")
        row = widgets.Box([msg, upload_button])

        def read_file(inputs):
            self.program = str(inputs["new"][0]["content"], "utf8")

        upload_button.observe(read_file, names="value")

        return row

    def viscoso_widget(self, sol):

        html = sol.log.to_viscoso_widget()

        class VisCoSoWidget(anywidget.AnyWidget):
            _esm = """
            function appendHtml(el, str) {
                var div = document.createElement('div');
                div.innerHTML = str;
                while (div.children.length > 0) {
                    el.appendChild(div.children[0]);
                    MathJax.Hub.Queue(["Typeset",MathJax.Hub,div.children[0]]);
                }
            }

            export function render(view) {
                let widget = () => view.model.get("value");
                appendHtml(view.el, widget());
            }
            """
            css_path = os.path.join(ROOT_DIR, "src", "VisCoSo", "css", "viscoso.css")
            _css = open(css_path, "r").read()
            value = traitlets.Unicode(html).tag(sync=True)

        vscw = VisCoSoWidget()
        display(vscw, display_id="viscoso")
