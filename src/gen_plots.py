import seaborn as sns
import os
import pandas
import matplotlib
import matplotlib.pyplot as plt

from util import ROOT_DIR

OUT_DIR = os.path.join(ROOT_DIR, "tests", "results")
BENCHMARKS = os.path.join(OUT_DIR, "bench_results.csv")
COSO_STATS = os.path.join(OUT_DIR, "subs_results.csv")
WIDTH = 430
PGF = False

pal = sns.color_palette("colorblind")
if PGF:
    matplotlib.use("pgf")
    matplotlib.rcParams.update(
        {
            "pgf.texsystem": "pdflatex",
            "font.family": "serif",
            "text.usetex": True,
            "pgf.rcfonts": False,
        }
    )


def autolabel(rects, fmt=".2f"):
    # attach some text labels
    for rect in rects:
        height = int(rect.get_height())
        rect.axes.annotate(
            f"{height}",
            xy=(rect.get_x() + rect.get_width() / 2.0, height),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center",
            va="bottom",
            zorder=10,
        )


def set_size(width_pt, fraction=1, subplots=(1, 1)):
    """Set figure dimensions to sit nicely in our document.

    Parameters
    ----------
    width_pt: float
            Document width in points
    fraction: float, optional
            Fraction of the width which you wish the figure to occupy
    subplots: array-like, optional
            The number of rows and columns of subplots.
    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    # Width of figure (in pts)
    fig_width_pt = width_pt * fraction
    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    golden_ratio = (5**0.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio * (subplots[0] / subplots[1])

    return (fig_width_in, fig_height_in)


########################
## Sat problems stats ##
########################


def plot_benchmarks_sat(df, fig, axis):

    sns.barplot(
        ax=axis,
        x="n_solutions",
        y="time",
        hue="solver",
        data=df,
        palette=pal,
        log=True,
    )

    axis.set_yscale("log", nonpositive="clip")
    axis.set(ylabel="Seconds")
    axis.set(xlabel="\\# solutions")

    axis.axhline(y=300, color="r", linestyle="--")

    ticks = [10**x for x in range(-2, 3)] + [300]
    axis.set_yticks(ticks)
    axis.set_yticklabels(ticks, va="center")
    for tick in axis.xaxis.get_major_ticks():
        tick.tick1line.set_visible(False)
        tick.tick2line.set_visible(False)
        tick.label1.set_visible(False)
        tick.label2.set_visible(False)
    axis.xaxis.tick_bottom()
    ticks = []
    for l in sorted(df["n_solutions"].unique()):
        if l <= 10000:
            ticks.append(str(l))
        else:
            ticks.append(("{:.2e}".format(l)))
    axis.set_xticklabels(ticks)
    axis.tick_params(axis="x", labelrotation=90)
    handles, labels = axis.get_legend_handles_labels()
    axis.legend(
        handles,
        ["CoSo", "ASP", "sharpSAT", "Essence"],
        bbox_to_anchor=(0.2, -0.5),
        loc="lower center",
        ncol=2,
        frameon=True,
        title="Framework",
        fontsize="x-small",
        title_fontsize="small",
        borderaxespad=0,
    )

    axis.set_title("Satisfiable problems")

    plt.subplots_adjust(bottom=0.35, left=0.1)

    if PGF:
        w, h = set_size(WIDTH)
        fig.set_size_inches(w, h)
        path = os.path.join(OUT_DIR, "sat.pgf")
        plt.savefig(path, bbox_inches="tight")


##########################
## Unsat problems stats ##
##########################


def plot_benchmarks_unsat(df, fig, axis):

    print(df)

    sns.barplot(ax=axis, x="benchmark", y="time", hue="solver", data=df, palette=pal)

    axis.set_yscale("log", nonpositive="clip")
    axis.set(ylabel="Seconds")

    axis.axhline(y=300, color="r", linestyle="--")

    ticks = [10**x for x in range(-2, 3)] + [300]
    axis.set_yticks(ticks)
    axis.set_yticklabels(ticks, va="center")
    for tick in axis.xaxis.get_major_ticks():
        tick.tick1line.set_visible(False)
        tick.tick2line.set_visible(False)
        tick.label1.set_visible(False)
        tick.label2.set_visible(False)
    axis.xaxis.tick_bottom()
    # xlabs = [f"unsat_{i}" for i, x in enumerate(axis.get_xticklabels()) ]
    # xlabs = ["ms 20/15", "pm 20/15", "sq 15/10", "sq 20/15", "cp 15/5", "cp 20/10"]
    xlabs = df["benchmark"].unique()
    axis.set_xticklabels(xlabs, va="center")
    axis.tick_params(axis="x", which="major", pad=5)
    handles, labels = axis.get_legend_handles_labels()
    lgd = axis.legend(
        handles,
        ["CoSo", "ASP", "sharpSAT", "Essence"],
        bbox_to_anchor=(1.22, 0.5),
        loc="center right",
        ncol=1,
        frameon=True,
        title="Framework",
        fontsize="x-small",
        title_fontsize="small",
        borderaxespad=0,
    )
    # axis.legend([],[], frameon=False)
    axis.set_title("Unsatisfiable problems")

    plt.subplots_adjust(left=0.15)

    if PGF:
        w, h = set_size(WIDTH)
        fig.set_size_inches(w, h)
        path = os.path.join(OUT_DIR, "unsat.pgf")
        plt.savefig(path, bbox_extra_artists=(lgd,), bbox_inches="tight")


############################
## CoSo subproblems stats ##
############################


def plot_coso_stats(df, fig, axis):

    tab_count = df.pivot_table(
        index=["n_subproblems"],
        columns="Benchmark",
        values="count",
        aggfunc="sum",
        fill_value=0,
    )
    tab_time = df.pivot_table(
        index=["n_subproblems"], values="mean", aggfunc="mean", fill_value=0
    )
    ax5 = tab_count.plot(
        kind="bar", stacked=True, color=sns.color_palette("colorblind"), width=0.9
    )

    labels_real = []
    labels_synth = []
    for i in range(0, len(tab_count)):
        labels_real.append(tab_count.iloc[i, 0])
        labels_synth.append(tab_count.iloc[i, 1])
    labels = labels_real + labels_synth

    for i, patch in enumerate(ax5.patches):
        x, y = patch.get_xy()
        x += patch.get_width() / 2
        y += patch.get_height() / 2
        y = 1 if y == 0.5 else y
        l = "" if labels[i] == 0 else labels[i]
        ax5.annotate(l, (x, y), ha="center", va="center", c="black")

    ax6 = ax5.twinx()
    ax7 = ax5.twinx()
    ax6.plot(ax5.get_xticks(), tab_time["mean"], color=sns.color_palette("Paired")[1])
    ax5.set_yscale("symlog", linthresh=10)
    ax6.set_yscale("log", nonpositive="clip")
    ax7.set_yscale("log", nonpositive="clip")

    ax7.set_ylim(ax6.get_ylim())

    ax5.grid(False)
    ax6.grid(False)
    ax7.grid(False)
    ax5.yaxis.set_ticks([], minor=True)
    ax6.yaxis.set_ticks([], minor=True)
    ax7.yaxis.set_ticks([], minor=True)
    ticks = [10**x for x in range(0, 3)]
    ax5.set_yticks(ticks)
    ax5.set_yticklabels(ticks, va="center")
    ticks = [10**x for x in range(-2, 3)]
    ax6.set_yticks(ticks)
    ax6.set_yticklabels(ticks, va="center")
    ax7.set(yticks=ax6.get_yticks(), yticklabels=[], ylim=ax6.get_ylim())

    for y in ax6.get_yticks():
        ax7.axhline(y=y, color=ax6.get_lines()[-1].get_c(), alpha=0.3, linestyle="-")

    ax5.set_xlabel("\\# subproblems")
    ax5.set_ylabel("\\# benchmarks")
    ax6.set_ylabel("Seconds")
    ax5.set_title("CoSo avg. runtime vs. \\#subproblems")

    ax5.set_zorder(ax7.get_zorder() + 1)
    ax6.set_zorder(ax7.get_zorder() + 1)
    ax5.set_frame_on(False)
    ax6.set_frame_on(False)

    ax7.spines["bottom"].set_color("0.5")
    ax7.spines["top"].set_color("0.5")
    ax7.spines["right"].set_color("0.5")
    ax7.spines["left"].set_color("0.5")

    plt.subplots_adjust(bottom=0.15, left=0.1)

    if PGF:
        w, h = set_size(WIDTH)
        fig.set_size_inches(w, h)
        path = os.path.join(OUT_DIR, "subproblems.pgf")
        plt.savefig(path, bbox_inches="tight")


def load_data():

    df_bench = pandas.read_csv(BENCHMARKS, sep="\t")
    df_coso = pandas.read_csv(COSO_STATS, sep="\t")

    for r1 in df_bench.index:
        for r2 in df_bench.index:
            if df_bench["benchmark"][r1] == df_bench["benchmark"][r2]:
                if df_bench["n_subproblems"][r1] == -1:
                    if df_bench["n_subproblems"][r2] != -1:
                        df_bench["n_subproblems"][r1] = df_bench["n_subproblems"][r2]
                if df_bench["n_solutions"][r1] == -1:
                    if df_bench["n_solutions"][r2] != -1:
                        df_bench["n_solutions"][r1] = df_bench["n_solutions"][r2]

    df_sat = df_bench[df_bench["n_solutions"] > 0]
    df_unsat = df_bench[df_bench["n_solutions"] == 0]

    return (df_sat, df_unsat, df_coso)


def plot(pgf=False):
    PGF = pgf
    sns.set_theme(style="darkgrid", color_codes=True)
    fig, axes = plt.subplots(nrows=1, ncols=3)
    df_sat, df_unsat, df_coso = load_data()
    plot_benchmarks_sat(df_sat, fig, axes[0])
    plot_benchmarks_unsat(df_unsat, fig, axes[1])
    plot_coso_stats(df_coso, fig, axes[2])
    if not PGF:
        plt.show()


if __name__ == "__main__":
    plot()
