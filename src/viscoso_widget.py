import anywidget
import traitlets
import ipywidgets as widgets
from IPython.display import display, clear_output

from .cola_parser import Parser
from .problem import EmptyException
from .util import *


def show_widget():
    Widget()


class Widget(object):
    def __init__(self):

        self.program = "A CoLa program"
        self.text_area = None
        self.current_output = None

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
            value=self.program,
            placeholder="Type a CoLa program",
            description="",
            disabled=False,
            layout=widgets.Layout( height='200px', min_height='100px', width='50%')
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
            self.program = self.text_area.value
            if self.current_output is not None:
                self.current_output.close()
                # with self.current_output:
                clear_output()
                self.input_layout()
                display(solve_button)

            parser = Parser(cola=self.program)
            try:
                parser.parse()
                sol = parser.problem.solve(debug=False)
                self.viscoso_widget(sol)
                self.current_output = output
                display(output)
            except EmptyException:
                with output:
                    print("Could not find a problem :(")

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
            _css = open(CSS_VISCOSO, "r").read()
            value = traitlets.Unicode(html).tag(sync=True)

        vscw = VisCoSoWidget()
        display(vscw, display_id="viscoso")
