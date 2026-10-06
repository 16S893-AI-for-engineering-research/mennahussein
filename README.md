# Menna Hussein — Personal Website

A multi-page personal academic website featuring animated visuals.

## Pages
- `index.html` — Introduction / bio
- `research.html` — Research Interests (launch vehicle failures, liquid rocket engine failure propagation, mitigation, adaptive trajectory control)
- `project.html` — Featured project (agentic MPC configuration for launch vehicle failure mitigation)
- `figure.html` — Regenerated Figure 3: MILP reconstruction, original/computed comparison, modeling provenance and numerical checks
- `devlog.html` — Dev Log / build history, kept in sync with the site's git history

## Structure
- `style.css` — shared styling and animations (animated gradient background, starfield, scroll reveals, hover glow, orbiting element)
- `script.js` — starfield generation, scroll-reveal via IntersectionObserver, card glow follow, progress bar animation, active nav highlighting

## Usage
Open `index.html` in a browser, or serve locally:

```bash
python3 -m http.server 8000
```

Then visit `http://localhost:8000`.

## Regenerated figure

See [`scripts/README.md`](scripts/README.md) for the lossless-convexification
experiment, inferred inputs, and known discrepancies with the paper. Rebuild:

```bash
python3 -m pip install -r scripts/requirements-figure.txt
python3 scripts/regenerate_figure.py
python3 -m unittest discover -s tests -v
```

The website displays generated local SVGs and does not require a plotting CDN.
The original PDF is needed only for optional re-extraction of reference graphics.

## Notes
All bio/project content is currently placeholder text and should be updated with real details.
