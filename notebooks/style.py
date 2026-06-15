import plotly.graph_objects as go
import plotly.io as pio

MOCHA = {
    "base": "#1e1e2e",
    "surface": "#313244",
    "text": "#cdd6f4",
    "subtext": "#a6adc8",
    "colors": [
        "#cba6f7",  # mauve
        "#89b4fa",  # blue
        "#94e2d5",  # teal
        "#a6e3a1",  # green
        "#f38ba8",  # red
        "#fab387",  # peach
        "#f9e2af",  # yellow
        "#89dceb",  # sky
        "#b4befe",  # lavender
        "#eba0ac",  # maroon
    ]
}

pio.templates["catppuccin"] = go.layout.Template(
    layout=go.Layout(
        paper_bgcolor=MOCHA["base"],
        plot_bgcolor=MOCHA["base"],
        font=dict(color=MOCHA["text"]),
        colorway=MOCHA["colors"],
        xaxis=dict(gridcolor=MOCHA["surface"], linecolor=MOCHA["surface"]),
        yaxis=dict(gridcolor=MOCHA["surface"], linecolor=MOCHA["surface"]),
    )
)