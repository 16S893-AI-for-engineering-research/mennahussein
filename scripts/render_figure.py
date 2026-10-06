"""Static, offline-friendly scientific figures generated from solver output."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    'font.size': 10, 'axes.labelsize': 11, 'axes.titlesize': 12,
    'svg.fonttype': 'none', 'svg.hashsalt': 'lossless-figure-3',
})


def save(fig, output, name):
    fig.savefig(output / name, bbox_inches='tight', facecolor='white',
                metadata={'Date': None, 'Creator': 'scripts/regenerate_figure.py'})
    plt.close(fig)


def trajectory_axes():
    fig, ax = plt.subplots(figsize=(5.5, 5.5), layout='constrained')
    ax.set(xlim=(-40000, 60000), ylim=(-10000, 90000),
           xlabel='Downrange (m)', ylabel='Altitude (m)')
    ax.set_xticks(np.arange(-40000, 60001, 20000))
    ax.set_yticks(np.arange(-10000, 90001, 10000))
    ax.ticklabel_format(style='sci', axis='both', scilimits=(4, 4), useMathText=True)
    ax.tick_params(labelsize=13)
    ax.xaxis.label.set_size(15)
    ax.yaxis.label.set_size(15)
    ax.xaxis.get_offset_text().set_size(12)
    ax.yaxis.get_offset_text().set_size(12)
    ax.set_aspect('equal', adjustable='box')
    ax.grid(True, linestyle=':', color='#bbb', linewidth=.7)
    return fig, ax


def render(bundle, output):
    s, ref = bundle['solution'], bundle['reference']
    p = np.asarray(s['position_m'])
    rp = np.asarray(ref['figure3']['position_m'])
    for overlay in (False, True):
        fig, ax = trajectory_axes()
        ax.plot(*p[0], 'o', mfc='none', mec='black', ms=8, label='Initial position', zorder=5)
        ax.plot(p[:, 0], p[:, 1], '.-', color='blue', ms=5, lw=1.2,
                label='MILP reconstruction' if overlay else 'Trajectory')
        ax.plot(*p[-1], 'x', color='red', ms=9, mew=1.5, label='Final position', zorder=5)
        ax.axhline(0, color='black', lw=1.2, label='Surface')
        if overlay:
            ax.plot(rp[:, 0], rp[:, 1], '--', color='#c75b0c', lw=1.2,
                    label='Published Fig. 3 (reference)')
        ax.legend(loc='upper left', fontsize=10 if overlay else 12, framealpha=1, fancybox=False)
        save(fig, output, 'comparison.svg' if overlay else 'regenerated-figure-3.svg')

    t = np.asarray(s['time_s'])
    v, u = np.asarray(s['velocity_m_s']), np.asarray(s['control_m_s2'])
    fig, axes = plt.subplots(5, 1, figsize=(8, 10), sharex=True, layout='constrained')
    axes[0].plot(t, p[:, 0] / 1000, label='Downrange', color='blue')
    axes[0].plot(t, p[:, 1] / 1000, label='Altitude', color='red')
    axes[0].set_ylabel('Position (km)')
    axes[1].plot(t, v[:, 0], color='blue')
    axes[1].plot(t, v[:, 1], color='red')
    axes[1].set_ylabel('Velocity (m/s)')
    # ZOH: stairs start at t=0 and run exactly to tf, not midpoint-shifted.
    axes[2].stairs(u[:, 0], t, color='blue', baseline=None, label='Downrange thrust')
    axes[2].stairs(u[:, 1], t, color='red', baseline=None, label='Vertical thrust')
    axes[2].set_ylabel('Thrust (m/s²)')
    axes[3].stairs(s['thrust_magnitude_m_s2'], t, color='black', baseline=None,
                   label='Actual Euclidean norm')
    axes[3].stairs(s['sigma_m_s2'], t, color='#777', linestyle=':', baseline=None,
                   label='Polygon slack σ')
    axes[3].axhline(2, color='red', linestyle='--', lw=1, label='Nominal minimum / maximum')
    axes[3].axhline(10, color='red', linestyle='--', lw=1)
    axes[3].set(ylabel='Thrust norm (m/s²)', ylim=(0, 11))
    axes[4].stairs(s['gravity_y_m_s2'], t, color='blue', baseline=None)
    axes[4].set(ylabel='Gravity (m/s²)', xlabel='Time (s)', xlim=(0, s['tf_s']))
    for ax in axes:
        ax.grid(True, linestyle=':', alpha=.5)
    for ax in (axes[0], axes[2], axes[3]):
        ax.legend(fontsize=8, loc='best')
    save(fig, output, 'profiles.svg')

    fig, ax = plt.subplots(figsize=(8, 3.8), layout='constrained')
    scan = bundle['time_search']
    pairs = [(e['tf_s'], e['cost_m_s']) for e in scan['evaluations'] if e['cost_m_s'] is not None]
    times, costs = np.asarray(pairs).T
    ax.plot(times, costs, '.-', color='#235ca8', ms=3, label='Fixed-time MILP solves')
    best = scan['best']
    ax.plot(best['tf_s'], best['cost_m_s'], 'o', color='#198754',
            label=f"Best found: {best['tf_s']:.1f} s")
    ax.plot(s['tf_s'], s['cost_m_s'], 'x', color='#d04b20', ms=8, mew=2,
            label='Paper-reported time: 490 s')
    ax.set(xlabel='Final time (s)', ylabel='Discrete cost Δt ∑σ (m/s)',
           title='Final-time search — reconstructed model')
    ax.grid(True, linestyle=':', alpha=.5)
    ax.legend(fontsize=9)
    save(fig, output, 'time-search.svg')

    fig, ax = plt.subplots(figsize=(8, 4.6), layout='constrained')
    ax.plot(rp[:, 0] / 1000, rp[:, 1] / 1000, '--', color='black', label='Published Fig. 3')
    for key, label, color in [('solution', 'Figure-consistent reconstruction', '#235ca8'),
                              ('inner_polygon_solution', 'Only polygon changed to inner', '#d27715'),
                              ('written_parameters_solution', 'Inner + 15 intervals + zero initial velocity', '#8c4aaa')]:
        positions = np.asarray(bundle[key]['position_m']) / 1000
        ax.plot(positions[:, 0], positions[:, 1], '.-', color=color, label=label, ms=3)
    ax.axhline(0, color='black', lw=.8)
    ax.set(xlabel='Downrange (km)', ylabel='Altitude (km)',
           title='Sensitivity to the text–figure discrepancies (all at 490 s)')
    ax.grid(True, linestyle=':', alpha=.5)
    ax.legend(fontsize=8)
    save(fig, output, 'alternatives.svg')
