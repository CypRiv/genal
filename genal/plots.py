from __future__ import annotations

import warnings

import matplotlib
import numpy as np
import pandas as pd
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextToPath
from plotnine import (
    aes,
    after_scale,
    coord_cartesian,
    element_blank,
    element_line,
    element_text,
    expand_limits,
    guide_legend,
    geom_abline,
    geom_linerange,
    geom_point,
    geom_segment,
    geom_vline,
    ggplot,
    guides,
    labs,
    scale_color_manual,
    scale_fill_manual,
    scale_x_continuous,
    scale_y_continuous,
    scale_y_reverse,
    theme,
    theme_bw,
    theme_get,
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


def _validate_figure_size(figure_size):
    if figure_size is None:
        return None
    if not isinstance(figure_size, (tuple, list)) or len(figure_size) != 2:
        raise ValueError("figure_size must be a tuple like (width, height) in inches.")
    fig_width, fig_height = figure_size
    if not isinstance(fig_width, (int, float)) or not isinstance(fig_height, (int, float)):
        raise ValueError("figure_size must be a tuple of two numbers (width, height) in inches.")
    if fig_width <= 0 or fig_height <= 0:
        raise ValueError("figure_size values must be > 0.")
    return float(fig_width), float(fig_height)


# Text sizes are in points, so that they keep their printed size in figures of any size. In figures smaller than the
# compact size of a plot (the smallest size in which its text fits), text and marker sizes are reduced in proportion to
# the figure, to at least MIN_TEXT_SCALE times their size and MIN_TEXT_SIZE points (smaller text is hard to read, and
# macOS's Helvetica fails to render below about 7 pixels with some matplotlib versions) and, in figures narrower than
# NARROW_FIGURE_WIDTH inches, legends are placed below the plot. Titles that do not fit are wrapped.
COMPACT_FIGURE_SIZE = (5, 3.5)
MIN_TEXT_SCALE = 0.6
MIN_TEXT_SIZE = 6
NARROW_FIGURE_WIDTH = 4.5


def _text_scale(figure_size, compact_size=COMPACT_FIGURE_SIZE):
    """Text scale of a figure: 1 if it is not smaller than the compact size of the plot, in proportion otherwise."""
    scale = min(figure_size[0] / compact_size[0], figure_size[1] / compact_size[1])
    return 1.0 if scale >= 1 else max(MIN_TEXT_SCALE, scale)


def _text_size(size, scale):
    """Size (points) of a text of `size` points in a figure with the given text scale."""
    return size if scale == 1 else max(MIN_TEXT_SIZE, size * scale)


def _text_width(text, size):
    """Width (inches) of a one-line text of `size` points in the font of the plots (current plotnine theme)."""
    with matplotlib.rc_context(theme_get().rcParams):
        width, _, _ = TextToPath().get_text_width_height_descent(text, FontProperties(size=size), ismath=False)
    return width / 72


def _fit_title(title, length, size, sep=" "):
    """
    Title wrapped on several lines if it is longer than `length` inches with a font of `size` points: between the
    parts separated by `sep`, and between words in parts that are too long.
    """
    lines = []
    for part in title.split(sep):
        if lines and _text_width(lines[-1] + sep + part, size) <= length:
            lines[-1] += sep + part
        elif sep != " " and _text_width(part, size) > length:
            lines += _fit_title(part, length, size).split("\n")
        else:
            lines.append(part)
    return "\n".join(lines)


def _figure_size_theme(figure_size, scale):
    """
    Line spacing of wrapped titles and, in a figure smaller than the compact size of the plot (scale < 1), text and
    legend sizes in proportion to the figure, with legends below the plot in narrow figures.
    """
    if scale == 1:
        return theme(axis_title=element_text(linespacing=1.1), plot_subtitle=element_text(linespacing=1.1))
    elements = dict(
        axis_title=element_text(size=_text_size(12, scale), linespacing=1.1),
        axis_text=element_text(size=_text_size(10, scale)),
        legend_title=element_text(size=_text_size(11, scale)),
        legend_text=element_text(size=_text_size(8.8, scale)),
        plot_subtitle=element_text(size=_text_size(10, scale), linespacing=1.1),
        legend_key_size=15.84 * scale,
    )
    if figure_size[0] < NARROW_FIGURE_WIDTH:
        elements.update(legend_position="bottom", legend_box="vertical")
    return theme(**elements)


def _mrpresso_outlier_ids(mrpresso_outliers) -> set[str]:
    """
    Return the SNP IDs flagged by the MR-PRESSO outlier test (stored by :meth:`Geno.MRpresso`) as a set of strings.
    Raises an error if MR-PRESSO has not been run.
    """
    if mrpresso_outliers is None:
        raise ValueError(
            "Use_mrpresso_data is set to True but MRpresso results not found. Please run MRpresso first."
        )
    return {str(x) for x in mrpresso_outliers}


def _wald_only_methods(methods, res):
    """
    If the MR results only contain the Wald ratio (single instrument), return ["Wald"] so that it is plotted
    instead of the requested methods. Otherwise return the requested methods.
    """
    wald_name = MR_METHODS_NAMES["Wald"]
    if not res.empty and (res["method"] == wald_name).all():
        print("The MR results contain only the Wald ratio (single instrument): plotting the Wald ratio estimate.")
        return ["Wald"]
    return methods


def mr_loo_plot(
    mr_loo_results,
    snps_per_page=40,
    page=None,
    top_influential=True,
    filename=None,
    exposure_name=None,
    outcome_name=None,
    figure_size=None,
    *,
    methods=None,
    use_mrpresso_data=None,
    mr_data=None,
    mrpresso_outliers=None,
    mrpresso_subset_data=None,
    mr_loo_config=None,
):
    """
    Create a leave-one-out forest plot using results stored from :meth:`Geno.MR_loo`.
    """
    if mr_loo_results is None:
        raise ValueError("Run MR_loo() before MR_loo_plot().")

    if not isinstance(mr_loo_results, (list, tuple)) or len(mr_loo_results) < 9:
        raise ValueError("Invalid MR_loo results; re-run MR_loo() before plotting.")

    (
        loo_df,
        all_row,
        exp_name,
        stored_outcome_name,
        method_key,
        method_display_name,
        odds,
        _heterogeneity,
        used_mrpresso,
    ) = mr_loo_results[:9]
    stored_config = (
        mr_loo_results[9]
        if len(mr_loo_results) >= 10 and isinstance(mr_loo_results[9], dict)
        else {}
    )
    if mr_loo_config is None:
        mr_loo_config = stored_config
    elif not isinstance(mr_loo_config, dict):
        raise TypeError("mr_loo_config must be a dict or None")

    if not isinstance(snps_per_page, int) or snps_per_page < 5:
        raise ValueError("snps_per_page must be an integer >= 5.")
    if page is not None and (not isinstance(page, int) or page < 1):
        raise ValueError("page must be an integer >= 1.")
    if loo_df.empty:
        raise ValueError("MR_loo() returned no results; nothing to plot.")

    exp_label = exp_name if exposure_name is None else exposure_name
    outcome_label = stored_outcome_name if outcome_name is None else outcome_name
    figure_size_norm = _validate_figure_size(figure_size)

    use_mrpresso_effective = (
        used_mrpresso if use_mrpresso_data is None else bool(use_mrpresso_data)
    )
    outlier_ids: set[str] = set()
    if use_mrpresso_effective:
        outlier_ids = _mrpresso_outlier_ids(mrpresso_outliers)

    from .MR_tools import ci_critical_value

    plot_df = loo_df.copy().reset_index(drop=True)
    all_b = float(all_row.get("b", np.nan))
    all_se = float(all_row.get("se", np.nan))
    # 95% CI critical values: same distribution as the p-value of the method (t for MR Egger and the mode methods)
    loo_nsnp = plot_df["nSNP"] if "nSNP" in plot_df.columns else pd.Series(np.nan, index=plot_df.index)
    plot_df["crit"] = [ci_critical_value(method_display_name, n) for n in loo_nsnp]
    all_crit = ci_critical_value(method_display_name, all_row.get("nSNP", np.nan))
    plot_df["_orig_order"] = np.arange(len(plot_df))
    plot_df["influence"] = (plot_df["b"] - all_b).abs()
    plot_df["estimate_sort"] = np.exp(plot_df["b"]) if odds else plot_df["b"]
    if outlier_ids and "SNP" in plot_df.columns:
        plot_df["mrpresso_status"] = np.where(
            plot_df["SNP"].astype(str).isin(outlier_ids), "outlier", "non-outlier"
        )

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
            np.exp(all_b - all_crit * all_se)
            if not np.isnan(all_b) and not np.isnan(all_se)
            else np.nan
        )
        all_ci_upper = (
            np.exp(all_b + all_crit * all_se)
            if not np.isnan(all_b) and not np.isnan(all_se)
            else np.nan
        )
        null_line = 1
        x_label = "Odds ratio"
        ci_lower_all = np.exp(plot_df["b"] - plot_df["crit"] * plot_df["se"])
        ci_upper_all = np.exp(plot_df["b"] + plot_df["crit"] * plot_df["se"])
    else:
        all_est = all_b
        all_ci_lower = (
            all_b - all_crit * all_se
            if not np.isnan(all_b) and not np.isnan(all_se)
            else np.nan
        )
        all_ci_upper = (
            all_b + all_crit * all_se
            if not np.isnan(all_b) and not np.isnan(all_se)
            else np.nan
        )
        null_line = 0
        x_label = "Causal estimate"
        ci_lower_all = plot_df["b"] - plot_df["crit"] * plot_df["se"]
        ci_upper_all = plot_df["b"] + plot_df["crit"] * plot_df["se"]

    def _normalize_extra_methods(extra):
        if extra is None:
            return []
        if isinstance(extra, str):
            extra_list = [extra]
        elif isinstance(extra, (list, tuple)):
            extra_list = list(extra)
        else:
            raise TypeError("methods must be None, a string, or a list/tuple of strings")

        out = []
        seen = {method_key}
        for method in extra_list:
            if not isinstance(method, str):
                raise TypeError("methods must be None, a string, or a list/tuple of strings")
            if method not in MR_METHODS_NAMES.keys() or method in ("all", "Wald"):
                warnings.warn(
                    f"{method} is not an appropriate MR method. MR methods can be IVW, WM, Egger... Please refer to the documentation for more."
                )
                continue
            if method in seen:
                continue
            seen.add(method)
            out.append(method)
        return out

    extra_method_keys = _normalize_extra_methods(methods)

    def _compute_method_summary(method: str, subset_data):
        from .MR_tools import MR_func

        if mr_data is None:
            raise ValueError(
                "You need to run query_outcome() before calling MR_loo_plot with methods."
            )

        action = int(mr_loo_config.get("action", 2))
        eaf_threshold = float(mr_loo_config.get("eaf_threshold", 0.42))
        nboot = int(mr_loo_config.get("nboot", 1000))
        penk = int(mr_loo_config.get("penk", 20))
        phi = float(mr_loo_config.get("phi", 1))
        cpus = int(mr_loo_config.get("cpus", 1))
        seed = mr_loo_config.get("seed")

        res, _df_mr = MR_func(
            mr_data,
            methods=[method],
            action=action,
            eaf_threshold=eaf_threshold,
            nboot=nboot,
            penk=penk,
            phi=phi,
            name_exposure=exp_name,
            cpus=cpus,
            subset_data=subset_data,
            seed=seed,
        )

        method_name = MR_METHODS_NAMES[method]
        if isinstance(method_name, (list, tuple)):
            method_name = method_name[0]

        if res.empty:
            return None

        row = res[res["method"] == method_name]
        if row.shape[0] != 1:
            return None

        return {
            "SNP": str(method_name),
            "Estimate": float(row["b"].values[0]),
            "SE": float(row["se"].values[0]),
            "crit": ci_critical_value(method_name, row["nSNP"].values[0]),
        }

    # "MR-PRESSO corrected" summary row: only when MR_loo ran on all instruments (otherwise "All instruments"
    # is already the outlier-removed estimate)
    corrected_row = None
    if use_mrpresso_effective and mrpresso_subset_data is not None and not used_mrpresso:
        corrected = _compute_method_summary(method_key, mrpresso_subset_data)
        if corrected is not None:
            corrected["SNP"] = "MR-PRESSO corrected"
            corrected_row = corrected

    summary_rows = [{"SNP": "All instruments", "Estimate": all_b, "SE": all_se, "crit": all_crit}]

    for method in extra_method_keys:
        method_summary = _compute_method_summary(method, subset_data=None)
        if method_summary is None:
            method_name = MR_METHODS_NAMES[method]
            if isinstance(method_name, (list, tuple)):
                method_name = method_name[0]
            warnings.warn(
                f"The {method_name} ({method}) method could not be computed and will be excluded from the plot."
            )
            continue
        summary_rows.append(method_summary)

    if corrected_row is not None:
        summary_rows.append(corrected_row)

    summary_df = pd.DataFrame(summary_rows)
    if odds:
        summary_df["CI_lower"] = np.exp(summary_df["Estimate"] - summary_df["crit"] * summary_df["SE"])
        summary_df["CI_upper"] = np.exp(summary_df["Estimate"] + summary_df["crit"] * summary_df["SE"])
        summary_df["Estimate"] = np.exp(summary_df["Estimate"])
    else:
        summary_df["CI_lower"] = summary_df["Estimate"] - summary_df["crit"] * summary_df["SE"]
        summary_df["CI_upper"] = summary_df["Estimate"] + summary_df["crit"] * summary_df["SE"]
    summary_df = summary_df.drop(columns="crit")
    summary_df["row_type"] = "overall"

    x_scale_kwargs = {}
    finite_ci = pd.concat(
        [ci_lower_all, ci_upper_all, summary_df["CI_lower"], summary_df["CI_upper"]],
        ignore_index=True,
    ).to_numpy()
    finite_ci = finite_ci[np.isfinite(finite_ci)]
    if finite_ci.size:
        x_focus_min = float(np.min(finite_ci))
        x_focus_max = float(np.max(finite_ci))
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
            df_page["CI_lower"] = np.exp(df_page["b"] - df_page["crit"] * df_page["se"])
            df_page["CI_upper"] = np.exp(df_page["b"] + df_page["crit"] * df_page["se"])
        else:
            df_page["Estimate"] = df_page["b"]
            df_page["CI_lower"] = df_page["b"] - df_page["crit"] * df_page["se"]
            df_page["CI_upper"] = df_page["b"] + df_page["crit"] * df_page["se"]

        plot_cols = ["SNP", "Estimate", "CI_lower", "CI_upper"]
        if "mrpresso_status" in df_page.columns:
            plot_cols.append("mrpresso_status")
        df_plot = df_page.loc[:, plot_cols].copy()
        df_plot["row_type"] = "instrument"
        if "mrpresso_status" not in df_plot.columns:
            df_plot["mrpresso_status"] = np.nan

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
                        "mrpresso_status": [np.nan],
                    }
                ),
                summary_df.assign(mrpresso_status=np.nan),
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
        if figure_size_norm is None:
            fig_width, fig_height = 10, height
        else:
            fig_width, fig_height = figure_size_norm
        text_scale = _text_scale((fig_width, fig_height), (COMPACT_FIGURE_SIZE[0], min(height, 0.25 * (n_rows_total + 3))))
        row_label_width = max(_text_width(label, _text_size(10, text_scale)) for label in y_labels)

        df_instruments = df_plot[df_plot["row_type"] == "instrument"]
        outlier_present = (
            use_mrpresso_effective
            and (df_instruments["mrpresso_status"] == "outlier").any()
        )
        if outlier_present:
            df_outliers = df_instruments[df_instruments["mrpresso_status"] == "outlier"]
            df_non_outliers = df_instruments[df_instruments["mrpresso_status"] != "outlier"]
        else:
            df_outliers = df_instruments.iloc[0:0]
            df_non_outliers = df_instruments

        df_overall = df_plot[df_plot["row_type"] == "overall"]

        # Confidence intervals as plain lines (no caps); overall estimates as diamonds with the same line width
        plot = (
            ggplot(df_plot, aes(x="Estimate", y="y_pos"))
            + geom_vline(xintercept=null_line, linetype="dashed", color="#9AA0A6", size=0.6)
            + geom_segment(
                data=df_non_outliers,
                mapping=aes(x="CI_lower", xend="CI_upper", yend="y_pos"),
                color="#9AA0A6",
                size=0.45,
                alpha=0.9,
            )
            + geom_segment(
                data=df_outliers,
                mapping=aes(x="CI_lower", xend="CI_upper", yend="y_pos"),
                color="#D62728",
                size=0.45,
                alpha=0.9 if outlier_present else 0.0,
            )
            + (
                geom_point(
                    data=df_instruments,
                    mapping=aes(fill="mrpresso_status"),
                    shape="o",
                    color="#111111",
                    stroke=0.5,
                    size=2.4 * text_scale,
                    alpha=0.95,
                )
                if outlier_present
                else geom_point(
                    data=df_instruments,
                    color="#111111",
                    size=2.4 * text_scale,
                    alpha=0.95,
                )
            )
            + geom_segment(
                data=df_overall,
                mapping=aes(x="CI_lower", xend="CI_upper", yend="y_pos"),
                color="#111111",
                size=0.45,
            )
            + geom_point(
                data=df_overall,
                shape="D",
                color="#111111",
                size=2.6 * text_scale,
            )
            + scale_y_continuous(
                breaks=y_breaks,
                labels=y_labels,
                limits=y_limits,
                expand=(0, 0),
            )
            + labs(
                x=x_label,
                y="",
                subtitle=_fit_title(subtitle, fig_width - row_label_width - 0.2, _text_size(10, text_scale), sep=" | ") if subtitle else subtitle,
            )
            + scale_x_continuous(**x_scale_kwargs)
            + theme(
                figure_size=(fig_width, fig_height),
                axis_title=element_text(size=12),
                axis_text=element_text(size=10),
                plot_subtitle=element_text(size=10, color="#444", ha="left"),
                panel_grid_major_y=element_blank(),
                panel_grid_minor_y=element_blank(),
                panel_grid_minor_x=element_blank(),
                panel_grid_major_x=element_line(color="#E6E6E6", size=0.3),
                axis_ticks_major_y=element_blank(),
                axis_ticks_minor_y=element_blank(),
            )
            + _figure_size_theme((fig_width, fig_height), text_scale)
        )

        if (
            use_mrpresso_effective
            and (
                df_plot.loc[df_plot["row_type"] == "instrument", "mrpresso_status"]
                == "outlier"
            ).any()
        ):
            plot += scale_fill_manual(
                values={"non-outlier": "#111111", "outlier": "#D62728"},
                breaks=["outlier", "non-outlier"],
                name="MR-PRESSO",
            )
            plot += guides(fill=guide_legend(override_aes={"size": 3 * text_scale, "alpha": 1}))
        else:
            plot += theme(legend_position="none")

        if not np.isnan(all_est):
            plot += geom_vline(
                xintercept=all_est,
                linetype="solid",
                color="#111111",
                size=0.5,
                alpha=0.8,
            )

        if filename is not None:
            if multi_page:
                out_path = f"{filename}_page{p:02d}.png"
            elif len(pages) == 1:
                out_path = f"{filename}.png"
            else:
                out_path = f"{filename}_page{p:02d}.png"
            plot.save(out_path, dpi=500, width=fig_width, height=fig_height, verbose=False)

        plots.append(plot)

    if len(plots) == 1:
        return plots[0]
    return PlotFolio(plots)


def mr_plot(
    mr_results,
    methods,
    exposure_name=None,
    outcome_name=None,
    filename=None,
    figure_size=None,
    use_mrpresso_data=False,
    mrpresso_outliers=None,
):
    """
    Creates and returns a scatter plot of individual SNP effects with lines representing
    different Mendelian Randomization (MR) methods.
    """
    if mr_results is None:
        raise ValueError(
            "You need to run an MR analysis with the MR method before calling the MR_plot function."
        )

    df_mr = mr_results[1].copy()
    res = mr_results[0]
    exposure_name = mr_results[2] if not exposure_name else exposure_name
    exposure_name = (
        "Effect on the exposure" if not exposure_name else f"Effect on {exposure_name}"
    )
    outcome_name = mr_results[3] if not outcome_name else outcome_name
    outcome_name = (
        "Effect on the outcome" if not outcome_name else f"Effect on {outcome_name}"
    )
    figure_size_norm = _validate_figure_size(figure_size) or (10, 6)
    text_scale = _text_scale(figure_size_norm)
    legend_right = figure_size_norm[0] >= NARROW_FIGURE_WIDTH  # Otherwise below the plot, taking height from the y title

    # Orient instruments so all SNPs increase the exposure (common MR plotting convention).
    mask_flip = df_mr["BETA_e"] < 0
    if mask_flip.any():
        df_mr.loc[mask_flip, ["BETA_e", "BETA_o"]] *= -1

    mrpresso_status = None
    if use_mrpresso_data:
        outlier_ids = _mrpresso_outlier_ids(mrpresso_outliers)
        if outlier_ids:
            if "SNP" in df_mr.columns:
                snp_ids = df_mr["SNP"].astype(str)
            else:
                snp_ids = df_mr.index.astype(str)
            is_outlier = snp_ids.isin(outlier_ids)
            if is_outlier.any():
                mrpresso_status = np.where(is_outlier, "outlier", "non-outlier")
                df_mr = df_mr.copy()
                df_mr["mrpresso_status"] = mrpresso_status

    # Error bars without caps (+/- 1 SE on each axis)
    plot = (
        ggplot(df_mr, aes("BETA_e", "BETA_o"))
        + geom_segment(
            aes(x="BETA_e-SE_e", xend="BETA_e+SE_e", yend="BETA_o"),
            color="gray",
            size=0.1,
        )
        + geom_linerange(
            aes(ymin="BETA_o-SE_o", ymax="BETA_o+SE_o"),
            color="gray",
            size=0.1,
        )
        + geom_abline(slope=0, intercept=0, color="black")
        + labs(
            x=_fit_title(exposure_name, 0.8 * figure_size_norm[0], _text_size(12, text_scale)),
            y=_fit_title(outcome_name, (0.85 if legend_right else 0.45) * figure_size_norm[1], _text_size(12, text_scale)),
        )
        + theme(
            axis_title=element_text(size=12),
            axis_text=element_text(size=10),
            figure_size=figure_size_norm,
        )
        + expand_limits(x=0)
    )

    if mrpresso_status is not None:
        non_outliers = df_mr[df_mr["mrpresso_status"] == "non-outlier"]
        outliers = df_mr[df_mr["mrpresso_status"] == "outlier"]
        # Outline of the fill color, also in the MR-PRESSO legend keys (which both layers draw); not in the methods legend
        plot += geom_point(
            data=non_outliers,
            mapping=aes(fill="mrpresso_status", color=after_scale("fill")),
            shape="o",
            size=0.2,
            stroke=0.5,
            show_legend={"color": False},
        )
        plot += geom_point(
            data=outliers,
            mapping=aes(fill="mrpresso_status", color=after_scale("fill")),
            shape="o",
            size=1.8,
            stroke=0.4,
            show_legend={"color": False},
        )
    else:
        plot += geom_point(color="black", size=0.2)

    if mrpresso_status is not None:
        plot += scale_fill_manual(
            values={"non-outlier": "#111111", "outlier": "#D62728"},
            breaks=["outlier", "non-outlier"],
            name="MR-PRESSO",
        )
        plot += guides(fill=guide_legend(override_aes={"size": 4 * text_scale, "alpha": 1}))

    lines = []
    for method in _wald_only_methods(methods, res):
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
        plot += labs(color="MR Methods")  # Otherwise named after the first color mapping (the outline of the points)
        if not legend_right:
            plot += guides(color=guide_legend(ncol=2, title=""))  # Method names need no legend title
    plot += _figure_size_theme(figure_size_norm, text_scale)

    if filename:
        plot.save(
            f"{filename}.png",
            dpi=500,
            width=figure_size_norm[0],
            height=figure_size_norm[1],
            verbose=False,
        )

    return plot


def mr_funnel(
    mr_results,
    methods="IVW",
    exposure_name=None,
    outcome_name=None,
    se_method="nome",
    symmetric_x=True,
    add_null_line=False,
    figure_size=(10, 6),
    filename=None,
    use_mrpresso_data=False,
    mrpresso_outliers=None,
):
    """
    Creates and returns a funnel plot of single-SNP causal estimates (Wald ratios).

    The plot displays the Wald ratio estimate (x-axis) versus its standard error (y-axis),
    with the y-axis reversed so more precise SNPs appear at the top.
    """
    if mr_results is None:
        raise ValueError(
            "You need to run an MR analysis with the MR method before calling the MR_funnel function."
        )

    df_mr = mr_results[1].copy()
    res = mr_results[0]
    _ = mr_results[2] if not exposure_name else exposure_name
    _ = mr_results[3] if not outcome_name else outcome_name
    figure_size_norm = _validate_figure_size(figure_size) or (10, 6)
    text_scale = _text_scale(figure_size_norm)

    # Orient instruments so all SNPs increase the exposure (common MR plotting convention).
    mask_flip = df_mr["BETA_e"] < 0
    if mask_flip.any():
        df_mr.loc[mask_flip, ["BETA_e", "BETA_o"]] *= -1

    beta_e = df_mr["BETA_e"].astype(float)
    beta_o = df_mr["BETA_o"].astype(float)
    se_e = df_mr["SE_e"].astype(float)
    se_o = df_mr["SE_o"].astype(float)

    ratio = beta_o / beta_e
    valid = (beta_e != 0) & np.isfinite(ratio)
    if not valid.any():
        raise ValueError("No valid single-SNP ratio estimates to plot.")

    if se_method == "nome":
        se_ratio = se_o / beta_e.abs()
    elif se_method == "delta":
        se_ratio = np.sqrt((se_o**2) / (beta_e**2) + (beta_o**2) * (se_e**2) / (beta_e**4))
    else:
        raise ValueError("se_method must be 'nome' or 'delta'.")

    valid = valid & np.isfinite(se_ratio) & (se_ratio > 0)
    if not valid.any():
        raise ValueError("No valid single-SNP ratio standard errors to plot.")

    mrpresso_status = None
    if use_mrpresso_data:
        outlier_ids = _mrpresso_outlier_ids(mrpresso_outliers)
        if outlier_ids:
            if "SNP" in df_mr.columns:
                snp_ids = df_mr["SNP"].astype(str)
            else:
                snp_ids = df_mr.index.astype(str)
            is_outlier = snp_ids.isin(outlier_ids)
            is_outlier_valid = is_outlier[valid]
            if is_outlier_valid.any():
                mrpresso_status = np.where(is_outlier_valid, "outlier", "non-outlier")

    plot_df = pd.DataFrame({"ratio": ratio[valid], "se_ratio": se_ratio[valid]}).reset_index(
        drop=True
    )
    if mrpresso_status is not None:
        plot_df["mrpresso_status"] = mrpresso_status

    # Normalize methods input: None -> no vertical lines; str -> single method; list/tuple -> multiple.
    if methods is None:
        method_keys = []
    elif isinstance(methods, str):
        method_keys = [methods]
    elif isinstance(methods, (list, tuple)):
        method_keys = list(methods)
    else:
        raise TypeError("methods must be None, a string, or a list/tuple of strings")
    if method_keys:
        method_keys = _wald_only_methods(method_keys, res)

    vlines = []
    seen_methods = set()
    for method in method_keys:
        if not isinstance(method, str):
            raise TypeError("methods must be None, a string, or a list/tuple of strings")
        if method not in MR_METHODS_NAMES.keys():
            warnings.warn(
                f"{method} is not an appropriate MR method. MR methods can be IVW, WM, Egger... Please refer to the documentation for more."
            )
            continue

        method_name = MR_METHODS_NAMES[method]
        if isinstance(method_name, (list, tuple)):
            method_name = method_name[0]

        if method_name in seen_methods:
            continue
        seen_methods.add(method_name)

        res_row = res[res.method == method_name]
        if res_row.shape[0] == 0:
            warnings.warn(
                f"The {method_name} ({method}) method was not included in the MR method call and will be excluded from the plot."
            )
            continue

        vlines.append(
            {
                "xintercept": float(res_row["b"].values[0]),
                "MR Methods": method_name,
            }
        )

    vline_df = pd.DataFrame(vlines)

    plot = (
        ggplot(plot_df, aes("ratio", "se_ratio"))
        + scale_y_reverse()
        + labs(
            x="Single-SNP estimate",
            y="Standard error",
        )
        + theme(
            axis_title=element_text(size=12),
            axis_text=element_text(size=10),
            figure_size=figure_size_norm,
        )
    )

    if mrpresso_status is not None:
        non_outliers = plot_df[plot_df["mrpresso_status"] == "non-outlier"]
        outliers = plot_df[plot_df["mrpresso_status"] == "outlier"]
        # Outline of the fill color, also in the MR-PRESSO legend keys (which both layers draw); not in the methods legend
        plot += geom_point(
            data=non_outliers,
            mapping=aes(fill="mrpresso_status", color=after_scale("fill")),
            shape="o",
            size=0.4,
            stroke=0.5,
            show_legend={"color": False},
        )
        plot += geom_point(
            data=outliers,
            mapping=aes(fill="mrpresso_status", color=after_scale("fill")),
            shape="o",
            size=2,
            stroke=0.4,
            show_legend={"color": False},
        )
    else:
        plot += geom_point(color="black", size=0.4)

    if mrpresso_status is not None:
        plot += scale_fill_manual(
            values={"non-outlier": "#111111", "outlier": "#D62728"},
            breaks=["outlier", "non-outlier"],
            name="MR-PRESSO",
        )
        plot += guides(fill=guide_legend(override_aes={"size": 4 * text_scale, "alpha": 1}))

    if not vline_df.empty:
        plot += geom_vline(
            aes(xintercept="xintercept", color="MR Methods"),
            data=vline_df,
            linetype="solid",
            size=0.8,
            alpha=0.95,
        )
        plot += labs(color="MR Methods")  # Otherwise named after the first color mapping (the outline of the points)
        if figure_size_norm[0] < NARROW_FIGURE_WIDTH:
            plot += guides(color=guide_legend(ncol=2, title=""))  # Method names need no legend title
    plot += _figure_size_theme(figure_size_norm, text_scale)

    if add_null_line:
        plot += geom_vline(xintercept=0, linetype="solid", color="#888888", size=0.5, alpha=0.8)

    if symmetric_x:
        if not vline_df.empty:
            center = float(vline_df["xintercept"].iloc[0])
        else:
            center = float(plot_df["ratio"].median())

        max_dev = float(np.nanmax(np.abs(plot_df["ratio"].to_numpy(dtype=float) - center)))
        if not vline_df.empty:
            max_dev = max(
                max_dev,
                float(np.nanmax(np.abs(vline_df["xintercept"].to_numpy(dtype=float) - center))),
            )
        if add_null_line:
            max_dev = max(max_dev, abs(center))

        if not np.isfinite(max_dev) or max_dev == 0:
            scale = float(np.nanmax(np.abs(plot_df["ratio"].to_numpy(dtype=float))))
            scale = max(abs(center), scale, 1.0)
            max_dev = 0.1 * scale

        margin = 0.05 * max_dev
        span = max_dev + margin
        plot += coord_cartesian(xlim=(center - span, center + span))
    elif add_null_line:
        plot += expand_limits(x=0)

    if filename:
        plot.save(
            f"{filename}.png",
            dpi=500,
            width=figure_size_norm[0],
            height=figure_size_norm[1],
            verbose=False,
        )

    return plot


def mr_forest(
    mr_results,
    methods,
    exposure_name=None,
    outcome_name=None,
    odds=False,
    filename=None,
    figure_size=None,
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
    figure_size_norm = _validate_figure_size(figure_size)

    plot_data = []
    for method in _wald_only_methods(methods, res):
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
    if plot_df.empty:
        raise ValueError(
            "None of the requested methods are available in the MR results: nothing to plot."
        )
    # 95% CI critical values: same distribution as the p-value of each method (t for MR Egger and the mode methods)
    from .MR_tools import ci_critical_value

    crit = np.array([ci_critical_value(m, n) for m, n in zip(plot_df["Method"], plot_df["nSNP"])])
    if odds:
        plot_df["CI_lower"] = np.exp(plot_df["Estimate"] - crit * plot_df["SE"])
        plot_df["CI_upper"] = np.exp(plot_df["Estimate"] + crit * plot_df["SE"])
        plot_df["Estimate"] = np.exp(plot_df["Estimate"])
        x_label = f"Odds Ratio for {outcome_name} per unit increase in {exposure_name}"
        null_line = 1
    else:
        plot_df["CI_lower"] = plot_df["Estimate"] - crit * plot_df["SE"]
        plot_df["CI_upper"] = plot_df["Estimate"] + crit * plot_df["SE"]
        x_label = f"Effect on {outcome_name} per unit increase in {exposure_name}"
        null_line = 0

    plot_df["Method"] = pd.Categorical(
        plot_df["Method"], categories=plot_df["Method"].values[::-1], ordered=True
    )

    figure_size_norm = figure_size_norm or (7, max(3.0, 1.4 + 0.5 * len(plot_df)))  # Height depending on the number of methods
    text_scale = _text_scale(figure_size_norm, (COMPACT_FIGURE_SIZE[0], 3))
    label_width = max(_text_width(method, _text_size(10, text_scale)) for method in plot_df["Method"])
    # Confidence intervals as plain lines (no caps) over the null line; no horizontal grid lines
    plot = (
        ggplot(plot_df, aes(x="Estimate", y="Method"))
        + geom_vline(xintercept=null_line, linetype="dashed", color="#9AA0A6", size=0.6)
        + geom_segment(aes(x="CI_lower", xend="CI_upper", yend="Method"), color="#111111", size=0.6)
        + geom_point(size=3 * text_scale, color="#111111")
        + theme(
            axis_title=element_text(size=12),
            axis_text=element_text(size=10),
            figure_size=figure_size_norm,
            panel_grid_major_y=element_blank(),
            panel_grid_minor_y=element_blank(),
            axis_ticks_major_y=element_blank(),
        )
        + labs(x=_fit_title(x_label, figure_size_norm[0] - label_width - 0.2, _text_size(12, text_scale)), y="")
        + _figure_size_theme(figure_size_norm, text_scale)
    )

    if filename:
        plot.save(
            f"{filename}.png",
            dpi=500,
            width=figure_size_norm[0],
            height=figure_size_norm[1],
            verbose=False,
        )

    return plot
