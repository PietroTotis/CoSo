import seaborn as sns
import os
import pandas
import matplotlib
import matplotlib.pyplot as plt

plt.style.use("science")

from .util import ROOT_DIR

OUT_DIR = os.path.join(ROOT_DIR, "tests", "results")
BENCHMARKS = os.path.join(ROOT_DIR, "tests", "benchmarks", "synthetic")
COSO_STATS_BENCH = os.path.join(OUT_DIR, "subs_results_benchmarks.csv")
COSO_STATS_EX = os.path.join(OUT_DIR, "subs_results_examples.csv")
WIDTH = 430
PGF = False

# pal = sns.color_palette("colorblind")


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


def plot_benchmarks_sat(df, fig, axis, save):

    sns.barplot(
        ax=axis,
        x="n_solutions",
        y="time",
        hue="solver",
        data=df,
        # palette=pal,
        log=True,
    )
    axis.set_axisbelow(True)
    axis.grid(b=True, which="major", axis="y", linestyle="-")

    axis.set_yscale("log", nonpositive="clip")
    axis.set_ylabel("Seconds (log)", fontsize=10)
    axis.set_xlabel("\\# solutions", fontsize=10)

    axis.axhline(y=300, color="r", linestyle="--")

    ticks = [10**x for x in range(-2, 3)] + [300]
    axis.set_yticks(ticks)
    axis.set_yticklabels(ticks, va="center")
    axis.yaxis.set_ticks([], minor=True)
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
    axis.tick_params(
        axis="x",  # changes apply to the x-axis
        which="both",  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
    )
    handles, labels = axis.get_legend_handles_labels()
    # ph = [plt.plot([], marker="", ls="")[0]]
    # handles = ph + ph + handles
    axis.legend(
        handles,
        ["CoSo", "ASP", "sharpSAT", "Essence"],
        bbox_to_anchor=(0.5, 1.2),
        loc="upper center",
        ncol=4,
        frameon=True,
        title="Framework",
        fontsize=8,
        title_fontsize=10,
        borderaxespad=0,
    )

    axis.set_title("Satisfiable problems", fontsize="10", x=0.5, y=1.22)
    axis.tick_params(axis="both", labelsize=8)

    # plt.subplots_adjust(bottom=0.35, left=0.1)

    if save:
        w, h = set_size(WIDTH)
        fig.set_size_inches(w + 0.1, h)
        path = os.path.join(OUT_DIR, "sat.pgf")
        plt.savefig(path, bbox_inches="tight")
        plt.close(fig)


##########################
## Unsat problems stats ##
##########################


def plot_benchmarks_unsat(df, fig, axis, save):

    sns.barplot(
        ax=axis,
        x="benchmark",
        y="time",
        hue="solver",
        data=df,
        # palette=pal
    )
    axis.set_axisbelow(True)
    axis.grid(b=True, which="major", axis="y", linestyle="-")

    axis.set_yscale("log", nonpositive="clip")
    axis.set_ylabel("Seconds (log)", fontsize="10")
    axis.set_xlabel("Benchmark", fontsize="10")

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
    axis.tick_params(axis="both", labelsize=8)
    handles, labels = axis.get_legend_handles_labels()
    lgd = axis.legend(
        handles,
        ["CoSo", "ASP", "sharpSAT", "Essence"],
        bbox_to_anchor=(0.5, 1.2),
        loc="upper center",
        ncol=4,
        frameon=True,
        title="Framework",
        fontsize=8,
        title_fontsize=10,
        borderaxespad=0,
    )
    # axis.legend([],[], frameon=False)
    axis.set_title("Unsatisfiable problems", fontsize=10, x=0.5, y=1.22)
    axis.yaxis.set_ticks([], minor=True)
    axis.tick_params(
        axis="x",  # changes apply to the x-axis
        which="both",  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
    )

    plt.subplots_adjust(left=0.15)

    if save:
        w, h = set_size(WIDTH)
        fig.set_size_inches(w, h)
        path = os.path.join(OUT_DIR, "unsat.pgf")
        plt.savefig(path, bbox_extra_artists=(lgd,), bbox_inches="tight")
        plt.close(fig)


############################
## CoSo subproblems stats ##
############################


def plot_coso_stats(df, fig, axis, save):

    count = df.groupby(["n_subproblems", "origin"])["count"].sum().unstack(fill_value=0)
    time = df.groupby(["n_subproblems"], as_index=False)["time"].mean()

    color_cycle = axis._get_lines.prop_cycler
    color_line = next(color_cycle)["color"]
    color_real = next(color_cycle)["color"]
    color_synth = next(color_cycle)["color"]

    count.plot(
        kind="bar",
        stacked=True,
        # palette=[color_real, color_synth],
        color=sns.color_palette([color_real, color_synth]),
        width=0.9,
        ax=axis,
    ).legend(loc="upper center", bbox_to_anchor=(0.5, 1))
    plt.xlim(-1, None)

    n = len(count)
    for i, patch in enumerate(axis.patches):
        c = 0 if i < n else 1
        r = i if i < n else i - n
        x, y = patch.get_xy()
        x += patch.get_width() / 2
        y += patch.get_height() / 2
        l = "" if count.iloc[r, c] == 0 else count.iloc[r, c]
        axis.annotate(l, (x, y), ha="center", va="center", c="black", fontsize=10)

    axis.tick_params(
        axis="x",  # changes apply to the x-axis
        which="both",  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
    )
    ax6 = axis.twinx()
    ax7 = axis.twinx()
    axis.tick_params(axis="both", labelsize=8)
    ax6.tick_params(axis="both", labelsize=8)
    ax6.plot(
        axis.get_xticks(),
        time["time"],
        color=color_line
        # color=sns.color_palette("colorblind")[2]
    )
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

    axis.set_xlabel("\\# subproblems", fontsize=10)
    axis.set_ylabel("\\# benchmarks (log)", fontsize=10)
    ax6.set_ylabel("Seconds (log)", fontsize=10)
    axis.set_title("CoSo avg. runtime vs. \\#subproblems", fontsize=10)

    axis.set_zorder(ax7.get_zorder() + 1)
    ax6.set_zorder(ax7.get_zorder() + 1)
    axis.set_frame_on(False)
    ax6.set_frame_on(False)

    ax7.spines["bottom"].set_color("0.5")
    ax7.spines["top"].set_color("0.5")
    ax7.spines["right"].set_color("0.5")
    ax7.spines["left"].set_color("0.5")

    axis.legend(loc=9, prop={"size": 10})
    plt.subplots_adjust(bottom=0.15, left=0.1)

    if save:
        w, h = set_size(WIDTH)
        fig.set_size_inches(w, h)
        path = os.path.join(OUT_DIR, "subproblems.pgf")
        plt.savefig(path, bbox_inches="tight")
        plt.close(fig)


#####################
## Growing domains ##
#####################


def plot_growing_doms(df, fig, axis, name):

    for ind in df.index:
        s = df["benchmark"][ind]
        df["benchmark"][ind] = 0
        for i in range(1, 11):
            if f"_{i}" in s:
                df["benchmark"][ind] = i

    sns.lineplot(
        ax=axis,
        x="benchmark",
        y="time",
        hue="solver",
        data=df,
        marker=".",
        # palette=pal,
    )

    axis.set_axisbelow(True)
    axis.grid(b=True, which="major", axis="y", linestyle="-")
    axis.tick_params(
        axis="x",  # changes apply to the x-axis
        which="minor",  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
    )
    axis.yaxis.set_ticks([], minor=True)

    # axis.set_yscale("symlog", linthresh=10)

    axis.axhline(y=300, color="r", linestyle="--")
    axis.set_xlabel(name, fontsize=10)

    if name == "P3":
        axis.set_ylabel("Seconds", fontsize=10)
        ticks = [1, 10, 50, 100, 200, 300]
        axis.set_yticks(ticks)
        axis.set_yticklabels(ticks, va="center", fontsize=10)
    axis.tick_params(axis="both", labelsize=8)
    # if name == "P5":
    #     handles, labels = axis.get_legend_handles_labels()
    #     fig.legend(
    #         handles,
    #         ["CoSo", "ASP", "sharpSAT", "Essence"],
    #         bbox_to_anchor=(0.2, -0.5),
    #         loc="lower center",
    #         ncol=2,
    #         frameon=True,
    #         title="Framework",
    #         fontsize="x-small",
    #         title_fontsize="small",
    #         borderaxespad=0,
    #     )
    # else:
    axis.get_legend().remove()


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
    # sns.set_theme(style="whitegrid", color_codes=True, palette=pal)
    df_sat, df_unsat, df_coso, df_real, df_growing = load_data(bench_dir)

    if pgf:
        matplotlib.use("pgf")
        matplotlib.rcParams.update(
            {
                "pgf.texsystem": "pdflatex",
                "font.family": "serif",
                "text.usetex": True,
                "pgf.rcfonts": False,
            }
        )
    if pgf:
        fig1, ax1 = plt.subplots()
        plot_benchmarks_sat(df_sat, fig1, ax1, pgf)
        fig2, ax2 = plt.subplots()
        plot_benchmarks_unsat(df_unsat, fig2, ax2, pgf)
        fig3, ax3 = plt.subplots()
        plot_coso_stats(df_coso, fig3, ax3, pgf)
    else:
        fig, axes = plt.subplots(nrows=1, ncols=3)
        plot_benchmarks_sat(df_sat, fig, axes[0], pgf)
        plot_benchmarks_unsat(df_unsat, fig, axes[1], pgf)
        plot_coso_stats(df_coso, fig, axes[2], pgf)
    fig_g, axes_g = plt.subplots(nrows=1, ncols=3, sharey=True)
    plot_growing_doms(df_growing[1], fig_g, axes_g[0], "P3")
    plot_growing_doms(df_growing[2], fig_g, axes_g[1], "P4")
    plot_growing_doms(df_growing[3], fig_g, axes_g[2], "P5")
    handles, labels = axes_g[2].get_legend_handles_labels()
    # ph = [plt.plot([], marker="", ls="")[0]]
    # handles = ph + handles
    # labels = ["Framework:"] + labels
    fig_g.legend(
        handles,
        ["CoSo", "ASP", "sharpSAT", "Essence"],
        bbox_to_anchor=(0.5, 1.0),
        loc="upper center",
        ncol=4,
        frameon=True,
        title="Framework",
        fontsize=8,
        title_fontsize=10,
        borderaxespad=0,
    )
    fig_g.tight_layout(pad=0.5)
    if pgf:
        w, h = set_size(WIDTH)
        fig_g.set_size_inches(w, h)
        path = os.path.join(OUT_DIR, "growing_domains.pgf")
        plt.savefig(path, bbox_inches="tight")
    else:
        plt.show()


if __name__ == "__main__":
    plot("")
