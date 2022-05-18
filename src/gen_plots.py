import seaborn as sns
import os
import pandas
import matplotlib
import matplotlib.pyplot as plt

from util import ROOT_DIR

OUT_DIR = os.path.join(ROOT_DIR, "tests", "results")
BENCHMARKS = os.path.join(ROOT_DIR, "tests", "benchmarks", "synthetic")
COSO_STATS_BENCH = os.path.join(OUT_DIR, "subs_results_benchmarks.csv")
COSO_STATS_EX = os.path.join(OUT_DIR, "subs_results_examples.csv")
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
    xlabs = [
        "cp 15/5",
        "cp 20/10",
        "ms 20/15",
        "pm 20/15",
        "sq 15/10",
        "sq 20/15",
    ]
    # xlabs = df["benchmark"].unique()
    axis.set_xticklabels(xlabs, va="top")
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

    count = df.groupby(["n_subproblems", "origin"])["count"].sum().unstack(fill_value=0)
    time = df.groupby(["n_subproblems"], as_index=False)["time"].mean()

    count.plot(
        kind="bar",
        stacked=True,
        color=sns.color_palette("colorblind"),
        width=0.9,
        ax=axis,
    )

    # labels_real = []
    # labels_synth = []
    # for i in range(0, len(df)):
    #     labels_real.append(df.iloc[i, 0])
    #     labels_synth.append(df.iloc[i, 1])
    # labels = labels_real + labels_synth

    n = len(count)
    for i, patch in enumerate(axis.patches):
        c = 0 if i < n else 1
        r = i if i < n else i - n
        x, y = patch.get_xy()
        x += patch.get_width() / 2
        y += patch.get_height() / 2
        # y = 0.8 if y == 0.5 else y
        l = "" if count.iloc[r, c] == 0 else count.iloc[r, c]
        axis.annotate(l, (x, y), ha="center", va="center", c="black", fontsize=10)

    ax6 = axis.twinx()
    ax7 = axis.twinx()
    ax6.plot(axis.get_xticks(), time["time"], color=sns.color_palette("colorblind")[2])
    axis.set_yscale("symlog", linthresh=10)
    ax6.set_yscale("log", nonpositive="clip")
    ax7.set_yscale("log", nonpositive="clip")

    ax7.set_ylim(ax6.get_ylim())

    axis.grid(False)
    ax6.grid(False)
    ax7.grid(False)
    axis.yaxis.set_ticks([], minor=True)
    ax6.yaxis.set_ticks([], minor=True)
    ax7.yaxis.set_ticks([], minor=True)
    ticks = [10**x for x in range(0, 3)]
    axis.set_yticks(ticks)
    axis.set_yticklabels(ticks, va="center")
    ticks = [10**x for x in range(-1, 3)]
    ax6.set_yticks(ticks)
    ax6.set_yticklabels(ticks, va="center")
    ax7.set(yticks=ax6.get_yticks(), yticklabels=[], ylim=ax6.get_ylim())

    for y in ax6.get_yticks():
        ax7.axhline(y=y, color=ax6.get_lines()[-1].get_c(), alpha=0.3, linestyle="-")

    axis.set_xlabel("\\# subproblems")
    axis.set_ylabel("\\# benchmarks")
    ax6.set_ylabel("Seconds")
    axis.set_title("CoSo avg. runtime vs. \\#subproblems")

    axis.set_zorder(ax7.get_zorder() + 1)
    ax6.set_zorder(ax7.get_zorder() + 1)
    axis.set_frame_on(False)
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


#####################
## growing domains ##
#####################


def plot_growing_doms(df, fig, axis, name, plot_legend=False):

    for ind in df.index:
        s = df["benchmark"][ind]
        df["benchmark"][ind] = 0
        for i in range(1, 7):
            if f"_{i}" in s:
                df["benchmark"][ind] = i

    sns.lineplot(
        ax=axis,
        x="benchmark",
        y="time",
        hue="solver",
        data=df,
        # palette=pal,
    )

    # axis.set_yscale("symlog", linthresh=10)

    if plot_legend:
        handles, labels = axis.get_legend_handles_labels()
        fig.legend(
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

    axis.axhline(y=300, color="r", linestyle="--")
    axis.set_xlabel(name)
    axis.set_ylabel("Seconds (log)")

    ticks = [10**x for x in range(0, 3)]
    axis.set_yticks(ticks)
    axis.set_yticklabels(ticks, va="center")

    if PGF:
        w, h = set_size(WIDTH)
        fig.set_size_inches(w, h)
        path = os.path.join(OUT_DIR, "subproblems.pgf")
        plt.savefig(path, bbox_inches="tight")


def load_folder(folder):
    b_name = f"bench_results_{folder}.csv"
    b_path = os.path.join(OUT_DIR, b_name)
    s_name = f"subs_results_{folder}.csv"
    s_path = os.path.join(OUT_DIR, s_name)
    if os.path.exists(b_path):
        df_bench = pandas.read_csv(b_path, sep="\t")
        df_subs = pandas.read_csv(s_path, sep="\t")
        for r1 in df_bench.index:
            for r2 in df_bench.index:
                if df_bench["benchmark"][r1] == df_bench["benchmark"][r2]:
                    if df_bench["n_subproblems"][r1] == -1:
                        if df_bench["n_subproblems"][r2] != -1:
                            df_bench["n_subproblems"][r1] = df_bench["n_subproblems"][
                                r2
                            ]
                    if df_bench["n_solutions"][r1] == -1:
                        if df_bench["n_solutions"][r2] != -1:
                            df_bench["n_solutions"][r1] = df_bench["n_solutions"][r2]
        return (df_bench, df_subs)
    else:
        p = pandas.DataFrame()
        return (p, p)


def load_data(bench_dir):

    # examples_dir = os.path.join(bench_dir, "examples")
    df_bench_real, df_coso_real = load_folder("examples")
    df_coso_synth = pandas.DataFrame()
    df_bench_synth = pandas.DataFrame()

    for folder in os.listdir(bench_dir):
        if folder != "examples":
            df_bench_folder, df_subs_folder = load_folder(folder)
            df_bench_synth = pandas.concat([df_bench_synth, df_bench_folder])
            df_coso_synth = pandas.concat([df_coso_synth, df_subs_folder])

    df_sat = df_bench_synth[df_bench_synth["n_solutions"] > 0]
    df_unsat = df_bench_synth[df_bench_synth["n_solutions"] == 0]
    # df_unknown = df_bench_synth[df_bench_synth["n_solutions"] == -1]

    df_coso_synth["origin"] = "synthetic"
    df_coso_real["origin"] = "real"
    df_coso = pandas.concat([df_coso_synth, df_coso_real])
    df_coso = df_coso[df_coso["n_subproblems"] >= 0]

    df_bench_growing = {}
    df_bench_growing_all, _ = load_folder("growing_domains")
    df_bench_growing[1] = df_bench_growing_all[
        df_bench_growing_all["benchmark"].str.startswith("h730")
    ]
    df_bench_growing[2] = df_bench_growing_all[
        df_bench_growing_all["benchmark"].str.startswith("m722")
    ]
    df_bench_growing[3] = df_bench_growing_all[
        df_bench_growing_all["benchmark"].str.startswith("m617")
    ]

    return (df_sat, df_unsat, df_coso, df_bench_real, df_bench_growing)


def plot(bench_dir, pgf=False):
    PGF = pgf
    sns.set_theme(style="darkgrid", color_codes=True)
    fig, axes = plt.subplots(nrows=1, ncols=3)
    df_sat, df_unsat, df_coso, df_real, df_growing = load_data(bench_dir)
    plot_benchmarks_sat(df_sat, fig, axes[0])
    plot_benchmarks_unsat(df_unsat, fig, axes[1])
    plot_coso_stats(df_coso, fig, axes[2])
    fig_g, axes_g = plt.subplots(nrows=1, ncols=3)
    plot_growing_doms(df_growing[1], fig_g, axes_g[0], "h730", plot_legend=True)
    plot_growing_doms(df_growing[2], fig_g, axes_g[1], "m722")
    plot_growing_doms(df_growing[3], fig_g, axes_g[2], "m617")
    if not PGF:
        plt.show()


if __name__ == "__main__":
    plot("")
