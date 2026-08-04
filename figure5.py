import matplotlib.pyplot as plt
import numpy as np
from ase.db import connect

# Database files
db_files = {
    'Terrace': ('dB/Au_T.db', 'dB/Pd_T.db', 'dB/notT.db'),
    'Edge':    ('dB/Au_E.db', 'dB/Pd_E.db', 'dB/notE.db'),
    'Kink':    ('dB/Au_K.db', 'dB/Pd_K.db', 'dB/notK.db'),
    'Ad':      ('dB/Au_A.db', 'dB/Pd_A.db', 'dB/notA.db'),
}

# Electrolytes and standard potentials
electrolytes = {
    'b) H$_{2}$O':  {'Au': 1.52,  'Pd': 0.92},
    'c) HCl':  {'Au': 1.001, 'Pd': 0.559},
    'd) HBr':  {'Au': 0.859, 'Pd': 0.507},
    'e) HI':   {'Au': 0.557, 'Pd': 0.245},
}

# Plot & label settings
colors = {
    'Terrace': 'teal',
    'Edge': 'lightcoral',
    'Kink': 'mediumaquamarine',
    'Ad': 'thistle'
}
labels = {
    'Terrace': 'Terrace (CN = 9)',
    'Edge': 'Edge (CN = 7)',
    'Kink': 'Kink (CN = 6)',
    'Ad': 'Ad (CN = 4)'
}
w = 0.020
fig, axs = plt.subplots(2, 2, figsize=(12, 8))
axs = axs.flatten()
legend_handles = []

# Plotting
for idx, (electrolyte, E0s) in enumerate(electrolytes.items()):
    ax = axs[idx]
    all_potentials = []
    label_potentials = []

    # Lists for all Au and Pd potentials combined (all sites)
    all_au_potentials = []
    all_pd_potentials = []

    for site, (au_db, pd_db, ref_db) in db_files.items():
        # Au potentials
        potentials_au = []
        for row in connect(au_db).select():
            SlabId = row.SlabId
            for ref in connect(ref_db).select(SlabId=SlabId):
                energy = ref.energy - 2.56 - row.energy
                potential = E0s['Au'] + energy / 3
                potentials_au.append(potential)
        potentials_au = np.array(potentials_au)
        if potentials_au.size > 0:
            all_au_potentials.extend(potentials_au)
            all_potentials.append(potentials_au)
            label_potentials.append((site, potentials_au, 'Au'))
            n, bins, patches = ax.hist(
                potentials_au,
                bins=np.arange(min(potentials_au), max(potentials_au) + w, w),
                alpha=0.85,
                color=colors[site],
                label=labels[site],
                histtype='stepfilled',
                edgecolor='black',
                zorder=5
            )
            if idx == 0:
                legend_handles.append(patches[0])

        # Pd potentials
        potentials_pd = []
        for row in connect(pd_db).select():
            SlabId = row.SlabId
            for ref in connect(ref_db).select(SlabId=SlabId):
                energy = ref.energy - 3.22 - row.energy
                potential = E0s['Pd'] + energy / 2
                potentials_pd.append(potential)
        potentials_pd = np.array(potentials_pd)
        if potentials_pd.size > 0:
            all_pd_potentials.extend(potentials_pd)
            all_potentials.append(potentials_pd)
            label_potentials.append((site, potentials_pd, 'Pd'))
            ax.hist(
                potentials_pd,
                bins=np.arange(min(potentials_pd), max(potentials_pd) + w, w),
                alpha=0.85,
                color=colors[site],
                histtype='stepfilled',
                edgecolor='black',
                zorder=5
            )

    # Draw background spans for Au and Pd potentials (all sites combined)
    if all_au_potentials:
        au_min = min(all_au_potentials)
        au_max = max(all_au_potentials)
        ax.axvspan(au_min, au_max, color='gold', alpha=0.05, label='_nolegend_', zorder=0)

    if all_pd_potentials:
        pd_min = min(all_pd_potentials)
        pd_max = max(all_pd_potentials)
        ax.axvspan(pd_min, pd_max, color='lightskyblue', alpha=0.085, label='_nolegend_', zorder=0)

    # Set x-range dynamically
    if all_potentials:
        combined = np.concatenate(all_potentials)
        min_pot = combined.min()
        max_pot = combined.max()
        ax.set_xlim(min_pot - 0.1, max_pot + 0.1)

        # Vertical standard lines
        ax.axvline(E0s['Au'], color='red', linestyle='--', lw=1.5)
        ax.axvline(E0s['Pd'], color='blue', linestyle='--', lw=1.5)

        # Annotate minimum potential, text placed to the right of Pd vertical line
        min_group = min(label_potentials, key=lambda x: x[1].min())
        min_site, min_data, min_elem = min_group
        min_val = min_data.min()

        margin = 0.05
        text_x = E0s['Pd'] - margin  # shift left of Pd line
        y_max = ax.get_ylim()[1]
        text_y = y_max * 0.9
        arrow_y = y_max * 0.15

        ax.annotate(
            f'{min_val:.2f} V',
            xy=(min_val, arrow_y),
            xytext=(text_x, text_y),
            fontsize=9,
            ha='right',  # align text right since it's on the left side
            arrowprops=dict(arrowstyle='->', lw=1),
            bbox=dict(boxstyle="round,pad=0.3", fc="gold", alpha=0.5),
            zorder=10
        )

    # Panel title
    ax.text(0.01, 0.975, electrolyte, transform=ax.transAxes,
            fontsize=11, weight='bold', verticalalignment='top')

    ax.set_ylim(0, 40)

    # Add frequency label only to first and third panels (idx 0 and 2)
    if idx in [0, 2]:
        ax.set_ylabel('Frequency', fontsize=12)
        ax.set_yticks([])  # Hide y-axis tick marks and numbers
    else:
        ax.set_yticks([])  # Hide y-axis ticks for others

    ax.set_xlabel('$E^\\circ_\\mathrm{diss}$ [V vs SHE]', fontsize=11)

# Global legend
fig.legend(
    handles=legend_handles,
    labels=labels.values(),
    loc='upper center',
    ncol=4,
    fontsize=12,
    frameon=False
)

plt.subplots_adjust(left=0.16, right=0.80, top=0.75, bottom=0.10, wspace=0.07, hspace=0.22)

#plt.savefig('Figure5_.png', bbox_inches='tight', dpi=1200, transparent=True)

plt.show()

