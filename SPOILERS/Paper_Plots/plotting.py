"""

plotting.py

Author: Jordan Mirocha
Affiliation: McGill University
Created on: Tue 17 May 2022 11:13:07 EDT

Description:

"""

import os
import numpy as np
try:
    from .data import PATH
except:
    PATH = '.'
import matplotlib.pyplot as pl
import matplotlib.gridspec as gridspec

bbox = dict(facecolor='none', edgecolor='k', fc='w',
    boxstyle='round,pad=0.3', alpha=0.6, zorder=1000)

class PlotSummary(object):
    """
    Create a "punchline plot" a la Fig. 18 in https://arxiv.org/abs/2108.07282.

    .. note :: This is so far intended for plotting constraints on Ts or Ts/TR,
        but could be generalized to make plots involving the ionized fraction,
        parameter constraints, etc.

    """
    def __init__(self, constraints):
        self.constraints = constraints

    @property
    def num_constraints(self):
        return len(self.constraints)

    @property
    def ticks_Ts_over_Tr(self):
        """
        Can manually set preferred tick marks in Ts/Tr to display.
        """
        if not hasattr(self, '_ticks_Ts_over_Tr'):
            self._ticks_Ts_over_Tr = np.arange(0, 0.6, 0.1)
        return self._ticks_Ts_over_Tr

    @ticks_Ts_over_Tr.setter
    def ticks_Ts_over_Tr(self, value):
        self._ticks_Ts_over_Tr = np.array(value)

    def _get_redshifts(self, redshifts):
        # Figure out what redshift(s) we're plotting
        if (redshifts is None):
            _redshifts = []
            for i, con in enumerate(self.constraints):
                _redshifts.extend(con.redshifts)

            redshifts = np.sort(np.unique(_redshifts))
            print("No 'redshifts' supplied: will use all ({}).".format(redshifts))
        elif type(redshifts) not in [list, tuple, np.ndarray]:
            redshifts = np.array([redshifts])

        return redshifts

    def get_adiabatic_limit(self, z):
        """
        Return kinetic temperature in adiabatically-cooling IGM.

        .. note :: Pulls from a lookup table generated with CosmoRec.

        """

        if not hasattr(self, '_inits'):
            fn = 'inits_planck_TTTEEE_lowl_lowE_best.txt'
            zarr, xe, Tk = np.loadtxt('{}/{}'.format(PATH, fn), unpack=True)
            self._inits = {'z': zarr, 'xe': xe, 'Tk': Tk}

        return np.interp(z, self._inits['z'], self._inits['Tk'])

    def get_labels(self, param, redshift):
        if param == 'Ts':
            labels = r'$\overline{T}_S(z=%.1f) \ [\mathrm{K}]$' % redshift, \
                r'$\overline{T}_S / T_{{\rm radio}}$'
        else:
            raise NotImplemented('help')

        return labels

    def get_axes(self, param, redshifts=None, num=1, orientation='vertical',
        show_priors=True, panel_width_long_dim=3, panel_width_short_dim=3,
        show_adiabatic_limit=True, show_grid=True, hspace=None):
        """
        Setup the plot window and axes according to how many models we have,
        how many redshifts, etc.

        Parameters
        ----------
        param : str
            Parameter of interest, e.g., 'Ts'.

        """

        redshifts = self._get_redshifts(redshifts)

        # Figure out dimensions and size of plot window
        dim_long = panel_width_long_dim * len(self.constraints)
        dim_short = panel_width_short_dim * len(redshifts)
        figsize = (dim_short, dim_long) if orientation == 'vertical' \
            else (dim_long, dim_short)

        nrows = len(self.constraints) if orientation == 'vertical' \
            else len(redshifts)
        ncols = len(self.constraints) if orientation == 'horizontal' \
            else len(redshifts)

        fig = pl.figure(num=num, figsize=figsize)
        spec = gridspec.GridSpec(ncols=ncols, nrows=nrows, figure=fig, hspace=hspace)

        # Create individual axes.
        axes = []
        axes_r = []
        for i, con in enumerate(self.constraints):
            axes_by_z = []
            axes_by_z_r = []
            for j, z in enumerate(redshifts):

                # Main axis label and twin axis label
                labels = self.get_labels(param, z)

                is_base = i == len(self.constraints) - 1
                is_twin = i == 0

                # Could generalize this if we want some panels to be bigger
                # than others.
                if orientation == 'vertical':
                    axes_by_z.append(fig.add_subplot(spec[i,j]))
                    if show_grid:
                        axes_by_z[j].grid(color='k', linestyle='--',
                            linewidth=0.5, which='major', axis='x')
                else:
                    axes_by_z.append(fig.add_subplot(spec[j,i]))

                    if show_grid:
                        axes_by_z[j].grid(color='k', linestyle='--',
                            linewidth=0.5, which='major', axis='y')

                # Main axis: gets first label generated.
                if is_base:
                    if orientation == 'vertical':
                        axes_by_z[j].set_xlabel(labels[0])
                    else:
                        axes_by_z[j].set_ylabel(labels[0])
                else:
                    if orientation == 'vertical':
                        axes_by_z[j].set_xticklabels([])
                    else:
                        axes_by_z[j].set_yticklabels([])

                # Overplot
                if (param == 'Ts') and show_adiabatic_limit:
                    T = self.get_adiabatic_limit(z)
                    if orientation == 'vertical':
                        axes_by_z[j].fill_betweenx([-10, 10], 0, T, ls='-',
                            color='gray', facecolors='none', hatch='x')
                    else:
                        axes_by_z[j].fill_between([-10, 10], 0, T, ls='-',
                            color='gray', facecolors='none', hatch='x')

                # Add twin axis
                ax_r = axes_by_z[j].twinx() if orientation == 'horizontal' \
                    else axes_by_z[j].twiny()

                axes_by_z_r.append(ax_r)

            axes.append(axes_by_z)
            axes_r.append(axes_by_z_r)


        return fig, axes, axes_r

    def plot_limits(self, param, redshifts=None, num=1, fig=None, axes=None,
        orientation='vertical', show_priors=True, panel_width_long_dim=2,
        panel_width_short_dim=4, show_adiabatic_limit=True, two_tailed=True,
        conf_level='both', arrow_len=3, lw=1, capsize=1, twinax='top', annotate_kw={},
        plim=None, pticks=10, temp_ax='top', hspace=None):
        """


        """

        assert orientation == 'vertical', \
            "Only 'vertical' orientation possible right now."

        # Setup axes if user doesn't supply them
        if (fig is None) and (axes is None):
            fig, axes, axes_r = self.get_axes(param, redshifts=redshifts,
                num=num, orientation=orientation, show_priors=show_priors,
                panel_width_long_dim=panel_width_long_dim,
                panel_width_short_dim=panel_width_short_dim,
                show_adiabatic_limit=show_adiabatic_limit, hspace=hspace)

        redshifts = self._get_redshifts(redshifts)

        fmt = '|' if orientation == 'vertical' else '_'

        # Start plotting :)
        xlo = None
        xhi = None
        subvars = list('abcdefghijklmnop')
        for i, con in enumerate(self.constraints):

            if 'variant_ids' in con.data:
                pdf_labels = con.data['variant_ids']
            else:
                pdf_labels = None

            for j, z in enumerate(redshifts):

                ax = axes[i][j]
                ax_r = axes_r[i][j]

                # Annotate model name in corner
                if ('model_name' in con.data) and (j == len(redshifts) - 1):
                    ax.annotate(con.data['model_name'], (0.98, 0.9),
                        xycoords='axes fraction', ha='right', va='top',
                        bbox=bbox, **annotate_kw)

                # Start with main constraint
                idstr = ['1']
                is_var = [True]
                is_subvar = [False]
                all_lims = [con.get_limits(param, z, two_tailed=two_tailed)]
                for k, var in enumerate(con.variants):
                    if k == 0:
                        # Did just above
                        continue

                    all_lims.append(con.get_limits(param+'_{}'.format(var), z,
                        two_tailed=two_tailed))
                    idstr.append(var)

                    is_subvar.append(not var[-1].isdigit())
                    is_var.append(not is_subvar[-1])

                ct = -1
                ticks = []
                for k, lims in enumerate(all_lims):
                    # Compute 68 and 95% lower limits.
                    # Assumes 'lims' is (95% lower limit, 68% LL)

                    if is_subvar[k]:
                        dv = 0.25 * subvars.index(idstr[k][-1])
                    else:
                        dv = 0
                        ct += 1

                    ticks.append(ct + dv)

                    if 'colors' in con.data.keys():
                        colors = con.data['colors']
                    else:
                        colors = ['k'] * len(con.data['variants'])

                    # Plot the constraints
                    if orientation == 'vertical':
                        if conf_level == 68 and lims[1] is not None:
                            ax.errorbar(lims[1], ct+dv, xerr=arrow_len,
                                lolims=True, xlolims=True, fmt=fmt, lw=lw, capsize=capsize, markeredgewidth=lw,
                                color=colors[k])
                        elif conf_level == 95 and lims[0] is not None:
                            ax.errorbar(lims[0], ct+dv, xerr=arrow_len,
                                lolims=True, xlolims=True, fmt=fmt, lw=lw, capsize=capsize, markeredgewidth=lw,
                                color=colors[k])
                        elif conf_level == 'both':
                            if lims[0] is not None:
                                ax.errorbar(lims[0], ct+dv, xerr=arrow_len,
                                    lolims=True, xlolims=True, fmt=fmt, lw=lw, capsize=capsize, markeredgewidth=lw,
                                    color=colors[k])
                            if lims[1] is not None:
                                ax.scatter(lims[1], ct+dv, marker='4',
                                    color=colors[k], s=200, lw=lw)

                    else:
                        pass

                # Adjust limits of axes in parameter space direction
                # and determine ticks
                if plim is not None:
                    lo, hi = plim

                    if orientation == 'vertical':
                        ax.set_xlim(lo, hi)
                        ax_r.set_xlim(lo, hi)
                        ticks_r = self.ticks_Ts_over_Tr * 2.725 * (1 + z)
                        ticks_r_s = \
                            ['%.2f' % lab for lab in self.ticks_Ts_over_Tr]
                        ticks_t = np.arange(lo, hi+pticks, pticks)
                        ticks_t_s = \
                            ['%i' % lab for lab in ticks_t]
                    else:
                        ticks_r = ticks_r_s = []
                        ax.set_ylim(lo, hi)
                        ax_r.set_ylim(lo, hi)

                # Fix ticks along shared axes for different models
                if i == len(self.constraints) - 1:
                    # Bottom panel
                    if orientation == 'vertical':
                        if temp_ax == 'bottom':
                            pass
                        else:
                            xlim_here = ax.get_xlim()
                            ax.set_xticks(ticks_t, [])
                            ax.set_xlim(*xlim_here)
                            ax.set_xlabel('')
                            ax.xaxis.set_ticks_position('top')
                            xlim_here = ax_r.get_xlim()
                            ax_r.set_xticks(ticks_r, ticks_r_s)
                            ax_r.set_xlim(*xlim_here)
                            ax_r.set_xlabel(r'$\overline{T}_S / T_{{\rm radio}}$')
                            ax_r.xaxis.set_ticks_position('bottom')
                            ax_r.xaxis.set_label_position('bottom')
                    else:
                        ax_r.set_ylabel(r'$\overline{T}_S / T_{{\rm radio}}$')


                elif i == 0:
                    # This is the top panel
                    if orientation == 'vertical':
                        if temp_ax == 'bottom':
                            xlim_here = ax_r.get_xlim()
                            ax_r.set_xticks(ticks_r, ticks_r_s)
                            ax_r.set_xlim(*xlim_here)
                            ax_r.set_xlabel(r'$\overline{T}_S / T_{{\rm radio}}$')
                        else:
                            xlim_here = ax_r.get_xlim()
                            ax_r.set_xticks(ticks_r, [])
                            ax_r.set_xlim(*xlim_here)
                            ax_r.set_xlabel('')
                            ax_r.xaxis.set_ticks_position('bottom')
                            xlim_here = ax.get_xlim()
                            ax.set_xticks(ticks_t, ticks_t_s)
                            ax.set_xlim(*xlim_here)
                            ax.set_xlabel(r'$\overline{T}_S / \rm{K}$')
                            ax.xaxis.set_ticks_position('top')
                            ax.xaxis.set_label_position('top')
                    else:
                        ax_r.set_ylabel(r'$\overline{T}_S / T_{{\rm radio}}$')
                else:
                    if orientation == 'vertical':
                        xlim_here = ax.get_xlim()
                        ax.set_xticks(ticks_t, [])
                        ax.set_xlim(*xlim_here)
                        xlim_here = ax_r.get_xlim()
                        ax_r.set_xticks(ticks_r, [])
                        ax_r.set_xlim(*xlim_here)

                        if temp_ax == 'bottom':
                            pass
                        else:
                            ax.xaxis.set_ticks_position('top')
                            ax_r.xaxis.set_ticks_position('bottom')

                    else:
                        ax_r.set_yticklabels([])

                # Fix tick mark for shared redshift axes labels
                if j > 0:
                    if orientation == 'vertical':
                        ax.set_yticks(ticks, [])
                        #ax.set_yticklabels([])
                        ax.set_ylim(-1, sum(is_var))
                    else:
                        ax.set_xticks(ticks, [])
                        #ax.set_xticklabels([])
                        ax.set_xlim(-1, sum(is_var))

                # Replace tick mark labels with strings.
                if j == 0 and (orientation == 'vertical'):
                    ax.set_yticks(ticks, pdf_labels)
                    #ax.set_yticklabels(pdf_labels)
                    ax.set_ylim(-1, sum(is_var))
                else:
                    pass

        return fig, axes, axes_r
