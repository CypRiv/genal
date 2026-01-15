from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from plotnine import (
    aes,
    element_blank,
    element_line,
    element_text,
    expand_limits,
    geom_abline,
    geom_errorbar,
    geom_errorbarh,
    geom_point,
    geom_vline,
    ggplot,
    labs,
    scale_color_manual,
    scale_x_continuous,
    scale_y_continuous,
    theme,
    theme_bw,
)

from .constants import MR_METHODS_NAMES


class PlotFolio(list):
    """
    A list of plotnine ggplot objects that renders nicely in Jupyter Notebooks.
    """

    def _repr_html_(self) -> str:  # pragma: no cover
        import base64
        import io

        import matplotlib.pyplot as plt

        parts = [
            "<div style='display:flex; flex-direction:column; gap:18px'>",
            f"<div style='font-weight:600'>Leave-one-out MR plots ({len(self)} pages)</div>",
        ]
        for i, plot in enumerate(self, start=1):
            fig = plot.draw()
            buf = io.BytesIO()
            fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
            plt.close(fig)
            data = base64.b64encode(buf.getvalue()).decode("ascii")
            parts.append(
                f"<div><div style='color:#666; font-size:12px; margin:0 0 6px 0'>Page {i}</div>"
                f"<img style='max-width:100%; height:auto; border:1px solid #eee; border-radius:6px' "
                f"src='data:image/png;base64,{data}'/></div>"
            )
        parts.append("</div>")
        return "\n".join(parts)


def mr_loo_plot(
    mr_loo_results,
    snps_per_page=40,
    page=None,
    top_influential=True,
    filename=None,
    exposure_name=None,
    outcome_name=None,
):
    """
    Create a leave-one-out forest plot using results stored from :meth:`Geno.MR_loo`.
    """
    if mr_loo_results is None:
        raise ValueError("Run MR_loo() before MR_loo_plot().")

    (
        loo_df,
        all_row,
        exp_name,
        stored_outcome_name,
        _method_key,
        method_display_name,
        odds,
        _heterogeneity,
        used_mrpresso,
    ) = mr_loo_results

    if not isinstance(snps_per_page, int) or snps_per_page < 5:
        raise ValueError("snps_per_page must be an integer >= 5.")
    if page is not None and (not isinstance(page, int) or page < 1):
        raise ValueError("page must be an integer >= 1.")
    if loo_df.empty:
        raise ValueError("MR_loo() returned no results; nothing to plot.")

    exp_label = exp_name if exposure_name is None else exposure_name
    outcome_label = stored_outcome_name if outcome_name is None else outcome_name

    plot_df = loo_df.copy().reset_index(drop=True)
    all_b = float(all_row.get("b", np.nan))
    all_se = float(all_row.get("se", np.nan))
    plot_df["_orig_order"] = np.arange(len(plot_df))
    plot_df["influence"] = (plot_df["b"] - all_b).abs()
    plot_df["estimate_sort"] = np.exp(plot_df["b"]) if odds else plot_df["b"]

    if top_influential:
        plot_df = plot_df.sort_values(
            ["influence", "_orig_order"],
            ascending=[False, True],
            kind="mergesort",
        ).head(snps_per_page)
        plot_df = plot_df.sort_values(
            "estimate_sort", ascending=False, kind="mergesort", na_position="last"
        )
        pages = [1]
        total_pages = 1
    else:
        plot_df = plot_df.sort_values(
            "estimate_sort", ascending=False, kind="mergesort", na_position="last"
        )
        n_snps = len(plot_df)
        n_pages = int(np.ceil(n_snps / snps_per_page))
        total_pages = n_pages
        if page is not None:
            if page > n_pages:
                raise ValueError(f"Invalid page {page}; max page is {n_pages}.")
            pages = [page]
        else:
            pages = list(range(1, n_pages + 1))

    subtitle = (
        f"Method: {method_display_name} | Exposure: {exp_label} | Outcome: {outcome_label} | "
        f"N instruments: {len(loo_df)}"
    )
    if used_mrpresso:
        subtitle += " | MR-PRESSO outliers removed"

    if top_influential:
        print(
            f"MR_loo_plot: showing top {len(plot_df)} influential instruments (ordered by estimate). "
            f"Set top_influential=False to plot all {len(loo_df)} instruments with pagination."
        )
    elif total_pages > 1:
        if page is None:
            print(
                f"MR_loo_plot: returning {total_pages} pages (ordered by estimate). "
                f"Use page=1..{total_pages} to render a single page."
            )
            print("Tip: in Jupyter, the returned object renders all pages inline.")
        else:
            print(
                f"MR_loo_plot: rendering page {page}/{total_pages} (ordered by estimate). "
                f"Set page=None to render all pages."
            )

    if odds:
        all_est = np.exp(all_b) if not np.isnan(all_b) else np.nan
        all_ci_lower = (
            np.exp(all_b - 1.96 * all_se)
            if not np.isnan(all_b) and not np.isnan(all_se)
            else np.nan
        )
        all_ci_upper = (
            np.exp(all_b + 1.96 * all_se)
            if not np.isnan(all_b) and not np.isnan(all_se)
            else np.nan
        )
        null_line = 1
        x_label = "Odds ratio"
        ci_lower_all = np.exp(plot_df["b"] - 1.96 * plot_df["se"])
        ci_upper_all = np.exp(plot_df["b"] + 1.96 * plot_df["se"])
    else:
        all_est = all_b
        all_ci_lower = (
            all_b - 1.96 * all_se
            if not np.isnan(all_b) and not np.isnan(all_se)
            else np.nan
        )
        all_ci_upper = (
            all_b + 1.96 * all_se
            if not np.isnan(all_b) and not np.isnan(all_se)
            else np.nan
        )
        null_line = 0
        x_label = "Causal estimate"
        ci_lower_all = plot_df["b"] - 1.96 * plot_df["se"]
        ci_upper_all = plot_df["b"] + 1.96 * plot_df["se"]

    x_scale_kwargs = {}
    finite_ci = pd.concat([ci_lower_all, ci_upper_all], ignore_index=True).to_numpy()
    finite_ci = finite_ci[np.isfinite(finite_ci)]
    if finite_ci.size:
        x_focus_min = float(np.min(finite_ci))
        x_focus_max = float(np.max(finite_ci))
        if np.isfinite(all_ci_lower):
            x_focus_min = min(x_focus_min, float(all_ci_lower))
            x_focus_max = max(x_focus_max, float(all_ci_upper))

        xlim_low = min(x_focus_min, null_line)
        xlim_high = max(x_focus_max, null_line)
        span = xlim_high - xlim_low
        if span == 0:
            span = abs(xlim_high) if xlim_high else 1.0
        pad = 0.05 * span
        xlim_low -= pad
        xlim_high += pad
        x_scale_kwargs = {"limits": (xlim_low, xlim_high)}

    plots = []
    multi_page = (not top_influential) and (total_pages > 1)
    for p in pages:
        if top_influential:
            df_page = plot_df.copy()
        else:
            start = (p - 1) * snps_per_page
            end = p * snps_per_page
            df_page = plot_df.iloc[start:end].copy()

        if odds:
            df_page["Estimate"] = np.exp(df_page["b"])
            df_page["CI_lower"] = np.exp(df_page["b"] - 1.96 * df_page["se"])
            df_page["CI_upper"] = np.exp(df_page["b"] + 1.96 * df_page["se"])
        else:
            df_page["Estimate"] = df_page["b"]
            df_page["CI_lower"] = df_page["b"] - 1.96 * df_page["se"]
            df_page["CI_upper"] = df_page["b"] + 1.96 * df_page["se"]

        df_page["crosses_null"] = (
            df_page["CI_lower"].notna()
            & df_page["CI_upper"].notna()
            & (df_page["CI_lower"] <= null_line)
            & (df_page["CI_upper"] >= null_line)
        )

        df_plot = df_page.loc[
            :, ["SNP", "Estimate", "CI_lower", "CI_upper", "crosses_null"]
        ].copy()
        df_plot["row_type"] = "instrument"

        df_plot = pd.concat(
            [
                df_plot,
                pd.DataFrame(
                    {
                        "SNP": [" "],
                        "Estimate": [np.nan],
                        "CI_lower": [np.nan],
                        "CI_upper": [np.nan],
                        "row_type": ["separator"],
                        "crosses_null": [np.nan],
                    }
                ),
                pd.DataFrame(
                    {
                        "SNP": ["All instruments"],
                        "Estimate": [all_est],
                        "CI_lower": [all_ci_lower],
                        "CI_upper": [all_ci_upper],
                        "row_type": ["overall"],
                        "crosses_null": [np.nan],
                    }
                ),
            ],
            ignore_index=True,
        )

        df_plot["y_pos"] = np.arange(len(df_plot), dtype=float)[::-1]
        y_scale_data = df_plot.sort_values("y_pos", ascending=True, kind="mergesort")
        y_breaks = y_scale_data["y_pos"].tolist()
        y_labels = y_scale_data["SNP"].astype(str).tolist()
        y_limits = (-0.5, float(len(df_plot)) - 0.5)

        n_rows_total = len(df_plot)
        height = max(6, min(0.25 * (n_rows_total + 3), 20))

        caption = None
        caption_parts = []
        if multi_page:
            caption_parts.append(f"Page {p}/{total_pages}")
        elif top_influential:
            caption_parts.append("Top influential instruments")
        caption_parts.append("SNPs ordered by estimate")
        caption = " • ".join(caption_parts)

        plot = (
            ggplot(df_plot, aes(x="Estimate", y="y_pos"))
            + geom_errorbarh(
                data=df_plot[df_plot["row_type"] == "instrument"],
                mapping=aes(xmin="CI_lower", xmax="CI_upper"),
                height=0,
                color="#9AA0A6",
                size=0.45,
                alpha=0.9,
            )
            + geom_point(
                data=df_plot[df_plot["row_type"] == "instrument"],
                mapping=aes(color="crosses_null"),
                size=2.4,
                alpha=0.95,
            )
            + geom_errorbarh(
                data=df_plot[df_plot["row_type"] == "overall"],
                mapping=aes(xmin="CI_lower", xmax="CI_upper"),
                height=0,
                color="#111111",
                size=0.8,
            )
            + geom_point(
                data=df_plot[df_plot["row_type"] == "overall"],
                color="#111111",
                size=3.0,
            )
            + geom_vline(xintercept=null_line, linetype="dashed", color="#9AA0A6", size=0.6)
            + scale_y_continuous(
                breaks=y_breaks,
                labels=y_labels,
                limits=y_limits,
                expand=(0, 0),
            )
            + scale_color_manual(
                values={True: "#9AA0A6", False: "#2C7FB8"},
                na_value="#9AA0A6",
            )
            + labs(
                x=x_label,
                y="",
                title="Leave-one-out MR",
                subtitle=subtitle,
                caption=caption,
            )
            + theme_bw()
            + scale_x_continuous(**x_scale_kwargs)
            + theme(
                figure_size=(10, height),
                axis_title_x=element_text(size=12),
                axis_text_x=element_text(size=10),
                axis_text_y=element_text(size=10),
                plot_title=element_text(size=16, weight="bold", ha="left"),
                plot_subtitle=element_text(size=10, color="#444", ha="left"),
                plot_caption=element_text(size=9, color="#666", ha="left"),
                panel_grid_major_y=element_blank(),
                panel_grid_minor_y=element_blank(),
                panel_grid_minor_x=element_blank(),
                panel_grid_major_x=element_line(color="#E6E6E6", size=0.3),
                axis_ticks_major_y=element_blank(),
                    axis_ticks_minor_y=element_blank(),
                    legend_position="none",
                )
        )

        if not np.isnan(all_est):
            plot += geom_vline(
                xintercept=all_est,
                linetype="solid",
                color="#111111",
                size=0.9,
                alpha=0.85,
            )

        if filename is not None:
            if multi_page:
                out_path = f"{filename}_page{p:02d}.png"
            elif len(pages) == 1:
                out_path = f"{filename}.png"
            else:
                out_path = f"{filename}_page{p:02d}.png"
            plot.save(out_path, dpi=500, width=10, height=height, verbose=False)

        plots.append(plot)

    if len(plots) == 1:
        return plots[0]
    return PlotFolio(plots)


def mr_plot(mr_results, methods, exposure_name=None, outcome_name=None, filename=None):
    """
    Creates and returns a scatter plot of individual SNP effects with lines representing
    different Mendelian Randomization (MR) methods.
    """
    if mr_results is None:
        raise ValueError(
            "You need to run an MR analysis with the MR method before calling the MR_plot function."
        )

    df_mr = mr_results[1]
    res = mr_results[0]
    exposure_name = mr_results[2] if not exposure_name else exposure_name
    exposure_name = (
        "Effect on the exposure" if not exposure_name else f"Effect on {exposure_name}"
    )
    outcome_name = mr_results[3] if not outcome_name else outcome_name
    outcome_name = (
        "Effect on the outcome" if not outcome_name else f"Effect on {outcome_name}"
    )

    df_mr["BETA_e"], df_mr["BETA_o"] = np.where(
        df_mr["BETA_e"] < 0,
        (-df_mr["BETA_e"], -df_mr["BETA_o"]),
        (df_mr["BETA_e"], df_mr["BETA_o"]),
    )

    plot = (
        ggplot(df_mr, aes("BETA_e", "BETA_o"))
        + geom_errorbarh(
            aes(xmin="BETA_e-SE_e", xmax="BETA_e+SE_e"),
            height=0,
            color="gray",
            size=0.1,
        )
        + geom_errorbar(
            aes(ymin="BETA_o-SE_o", ymax="BETA_o+SE_o"),
            width=0,
            color="gray",
            size=0.1,
        )
        + geom_point(color="black", size=0.2)
        + geom_abline(slope=0, intercept=0, color="black")
        + labs(x=exposure_name, y=outcome_name)
        + theme(
            axis_title=element_text(size=12),
            axis_text=element_text(size=10),
            figure_size=(10, 6),
        )
        + expand_limits(x=0)
    )

    lines = []
    for method in methods:
        if method not in MR_METHODS_NAMES.keys():
            warnings.warn(
                f"{method} is not an appropriate MR method. MR methods can be IVW, WM, Egger... Please refer to the documentation for more."
            )
            continue

        if not method.startswith("Egger"):
            method_name = MR_METHODS_NAMES[method]
            res_row = res[res.method == method_name]
            if res_row.shape[0] == 0:
                warnings.warn(
                    f"The {method_name} ({method}) method was not included in the MR method call and will be excluded from the plot."
                )
            elif res_row.shape[0] == 1:
                lines.append(
                    {"slope": res_row["b"].values[0], "intercept": 0, "MR Methods": method_name}
                )
        else:
            method_name = MR_METHODS_NAMES[method][0]
            method_name_intercept = MR_METHODS_NAMES[method][1]
            res_row = res[res.method == method_name]
            res_row_intercept = res[res.method == method_name_intercept]
            if res_row.shape[0] == 0:
                warnings.warn(
                    f"The {method_name} ({method}) method was not included in the MR method call and will be excluded from the plot."
                )
            elif res_row.shape[0] == 1 and res_row_intercept.shape[0] == 1:
                lines.append(
                    {
                        "slope": res_row["b"].values[0],
                        "intercept": res_row_intercept["b"].values[0],
                        "MR Methods": method_name,
                    }
                )

    if lines:
        line_data = pd.DataFrame(lines)
        plot += geom_abline(
            aes(slope="slope", intercept="intercept", color="MR Methods"), data=line_data
        )

    if filename:
        plot.save(f"{filename}.png", dpi=500, width=10, height=6, verbose=False)

    return plot


def mr_forest(
    mr_results,
    methods,
    exposure_name=None,
    outcome_name=None,
    odds=False,
    filename=None,
):
    """
    Creates and returns a forest plot of MR results, with one row per method.
    """
    if mr_results is None:
        raise ValueError(
            "You need to run an MR analysis with the MR method before calling the MR_forest function."
        )

    res = mr_results[0]
    exposure_name = mr_results[2] if not exposure_name else exposure_name
    exposure_name = "Exposure" if not exposure_name else exposure_name
    outcome_name = mr_results[3] if not outcome_name else outcome_name
    outcome_name = "Outcome" if not outcome_name else outcome_name

    plot_data = []
    for method in methods:
        if method not in MR_METHODS_NAMES.keys():
            warnings.warn(
                f"{method} is not an appropriate MR method. MR methods can be IVW, WM, Egger... Please refer to the documentation for more."
            )
            continue

        method_name = MR_METHODS_NAMES[method]
        if isinstance(method_name, tuple):
            method_name = method_name[0]

        res_row = res[res.method == method_name]
        if res_row.shape[0] == 0:
            warnings.warn(
                f"The {method_name} ({method}) method was not included in the MR method call and will be excluded from the plot."
            )
        elif res_row.shape[0] == 1:
            plot_data.append(
                {
                    "Method": method_name,
                    "Estimate": res_row["b"].values[0],
                    "SE": res_row["se"].values[0],
                    "nSNP": res_row["nSNP"].values[0],
                }
            )

    plot_df = pd.DataFrame(plot_data)
    if odds:
        plot_df["CI_lower"] = np.exp(plot_df["Estimate"] - 1.96 * plot_df["SE"])
        plot_df["CI_upper"] = np.exp(plot_df["Estimate"] + 1.96 * plot_df["SE"])
        plot_df["Estimate"] = np.exp(plot_df["Estimate"])
        x_label = f"Odds Ratio for {outcome_name} per unit increase in {exposure_name}"
        null_line = 1
    else:
        plot_df["CI_lower"] = plot_df["Estimate"] - 1.96 * plot_df["SE"]
        plot_df["CI_upper"] = plot_df["Estimate"] + 1.96 * plot_df["SE"]
        x_label = f"Effect on {outcome_name} per unit increase in {exposure_name}"
        null_line = 0

    plot_df["Method"] = pd.Categorical(
        plot_df["Method"], categories=plot_df["Method"].values[::-1], ordered=True
    )

    plot = (
        ggplot(plot_df, aes(x="Estimate", y="Method"))
        + geom_point(size=3)
        + geom_errorbarh(aes(xmin="CI_lower", xmax="CI_upper"), height=0.2)
        + theme(
            axis_title=element_text(size=12),
            axis_text=element_text(size=10),
            figure_size=(10, 6),
        )
        + labs(x=x_label, y="")
        + geom_vline(xintercept=null_line, linetype="dashed", color="gray")
    )

    if filename:
        plot.save(f"{filename}.png", dpi=500, width=10, height=6, verbose=False)

    return plot
