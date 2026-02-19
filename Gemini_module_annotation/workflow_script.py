
import pandas as pd
import urllib.parse
import os

def generate_module_report(input_csv, output_html, metadata_dict):
    """
    Generates an interactive HTML Functional Atlas from a module gene list CSV.
    
    Args:
        input_csv (str): Path to CSV with columns ['module_name', 'module_size', 'module_genes_hgnc']
        output_html (str): Path to save the generated report.
        metadata_dict (dict): Dictionary mapping module IDs to functional annotations.
    """
    
    # 1. Load Data
    df = pd.read_csv(input_csv)
    
    # 2. STRING URL Generator helper
    def get_string_url(genes_list, size_cutoff=500):
        genes = [g.strip() for g in genes_list.split(',') if g.strip()]
        if len(genes) > size_cutoff: return None
        display_genes = genes[:150]
        params = {
            "identifiers": "\r".join(display_genes),
            "species": 9606,
            "add_white_nodes": 0,
            "network_flavor": "confidence",
            "hide_labels": 0,
            "network_type": "functional",
            "required_score": 400
        }
        return f"https://string-db.org/api/image/network?{urllib.parse.urlencode(params)}"

    # 3. Custom SVG Histogram Generator
    def generate_enhanced_hist(df):
        sizes = df['module_size'].tolist()
        names = [n.replace("module_", "") for n in df['module_name']]
        max_size = max(sizes)
        # Scaled to accommodate labels
        svg_w, svg_h = 1000, 240
        margin_l, margin_r, margin_t, margin_b = 60, 20, 20, 70
        chart_w = svg_w - margin_l - margin_r
        chart_h = svg_h - margin_t - margin_b
        
        svg = f'<svg width="100%" height="{svg_h}" viewBox="0 0 {svg_w} {svg_h}" xmlns="http://www.w3.org/2000/svg" style="background:white; border-radius:8px;">'
        
        # Y-axis Scale
        y_ticks = [0, 500, 1000, 1500]
        for yt in y_ticks:
            y_pos = margin_t + chart_h - (yt / 1500) * chart_h
            svg += f'<line x1="{margin_l-5}" y1="{y_pos}" x2="{margin_l}" y2="{y_pos}" stroke="#94a3b8" />'
            svg += f'<text x="{margin_l-10}" y="{y_pos+4}" text-anchor="end" font-size="10" fill="#64748b">{yt}</text>'
            svg += f'<line x1="{margin_l}" y1="{y_pos}" x2="{svg_w-margin_r}" y2="{y_pos}" stroke="#f1f5f9" stroke-dasharray="4" />'
        
        svg += f'<line x1="{margin_l}" y1="{margin_t}" x2="{margin_l}" y2="{margin_t+chart_h}" stroke="#94a3b8" />'
        svg += f'<text x="{margin_l-45}" y="{margin_t+chart_h/2}" text-anchor="middle" font-size="10" fill="#64748b" transform="rotate(-90, {margin_l-45}, {margin_t+chart_h/2})">Gene Count</text>'
        
        bar_w = chart_w / len(sizes)
        for i, s in enumerate(sizes):
            h = (s / 1500) * chart_h
            x = margin_l + i * bar_w
            svg += f'<rect x="{x+2}" y="{margin_t+chart_h-h}" width="{bar_w-4}" height="{h}" fill="#3b82f6" rx="2"><title>Module {names[i]}: {s} genes</title></rect>'
            svg += f'<text x="{x + bar_w/2}" y="{margin_t+chart_h+15}" text-anchor="middle" font-size="8" fill="#94a3b8" transform="rotate(45, {x + bar_w/2}, {margin_t+chart_h+15})">{names[i]}</text>'
        
        svg += f'<text x="{margin_l + chart_w/2}" y="{svg_h - 10}" text-anchor="middle" font-size="12" font-weight="bold" fill="#64748b">Module</text>'
        svg += '</svg>'
        return svg

    # 4. Component Assembly
    hist_svg = generate_enhanced_hist(df)
    
    nav_links = ""
    summary_rows = ""
    module_cards = ""
    
    unique_genes_count = len(df['module_genes_hgnc'].str.split(',').explode().unique())
    avg_size = int(df['module_size'].mean())

    for _, row in df.iterrows():
        m_id = row['module_name']
        meta = metadata_dict.get(m_id, {"title": m_id, "short_theme": "Co-expression", "summary_theme": "N/A", "description": "No description available.", "cell_types": "N/A", "pathways": [], "disease": "N/A"})
        
        # Nav
        nav_links += f'<a href="#section_{m_id}"><span style="color:#94a3b8; font-size:0.7rem; margin-right:5px;">{m_id.split("_")[1]}</span> {meta["short_theme"]}</a>'
        
        # Summary Row
        summary_rows += f"<tr><td><strong>{m_id}</strong></td><td>{meta['title']}</td><td>{row['module_size']}</td><td>{meta['summary_theme'][:60]}...</td></tr>"
        
        # Details Card
        s_url = get_string_url(row['module_genes_hgnc'])
        full_genes = row['module_genes_hgnc'].replace(",", ", ")
        
        module_cards += f"""
        <div class="card" id="section_{m_id}">
            <div class="card-header">
                <span class="badge">{m_id}</span>
                <span class="size-label">Size: {row['module_size']} genes</span>
            </div>
            <div class="card-body">
                <div class="text-side">
                    <h2 class="module-title">{meta['title']}</h2>
                    <p class="summary"><strong>Primary Theme:</strong> {meta['summary_theme']}</p>
                    <p class="desc">{meta['description']}</p>
                    <div class="metadata-grid">
                        <div class="m-item"><strong>Cell Types:</strong> {meta['cell_types']}</div>
                        <div class="m-item"><strong>Clinical Link:</strong> {meta['disease']}</div>
                    </div>
                    <h3>Key Pathways</h3>
                    <ul class="pathway-list">
                        {"".join(f"<li>{p}</li>" for p in meta.get('pathways', []))}
                    </ul>
                </div>
                <div class="viz-side">
                    {f'<img src="{s_url}" class="string-img" alt="Network">' if s_url else '<div class="no-viz">Network omitted (large module)</div>'}
                </div>
            </div>
            <div class="gene-footer">
                <details>
                    <summary style="cursor:pointer; color:#3b82f6; font-weight:600;">View Full Gene List ({row['module_size']} genes)</summary>
                    <div style="margin-top:10px; line-height:1.4;">{full_genes}</div>
                </details>
            </div>
        </div>
        """

    # 5. Final Template
    full_html = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <title>Unsupervised Module Annotation</title>
        <style>
            :root {{ --primary: #0f172a; --accent: #3b82f6; --bg: #f8fafc; --text: #1e293b; --card: #ffffff; }}
            body {{ font-family: 'Inter', system-ui, sans-serif; background: var(--bg); color: var(--text); margin: 0; display: flex; }}
            nav {{ width: 240px; background: var(--primary); height: 100vh; position: fixed; overflow-y: auto; color: white; padding: 20px 10px; }}
            nav h3 {{ font-size: 0.9rem; text-transform: uppercase; letter-spacing: 0.1em; color: #94a3b8; margin-bottom: 20px; padding-left: 10px; }}
            nav a {{ display: flex; align-items: center; color: #cbd5e1; text-decoration: none; padding: 8px 12px; border-radius: 6px; font-size: 0.75rem; margin-bottom: 2px; transition: 0.2s; }}
            nav a:hover {{ background: #1e293b; color: white; }}
            main {{ flex: 1; margin-left: 260px; padding: 40px; max-width: 1100px; }}
            .dashboard-header {{ background: white; border: 1px solid #e2e8f0; border-radius: 12px; padding: 30px; margin-bottom: 40px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }}
            .subtitle {{ font-size: 1rem; color: #64748b; margin-top: -5px; margin-bottom: 25px; }}
            .stats {{ display: flex; gap: 40px; margin-bottom: 10px; }}
            .stat-box {{ flex: 1; }}
            .stat-val {{ font-size: 1.8rem; font-weight: 800; color: var(--accent); }}
            .stat-label {{ text-transform: uppercase; font-size: 0.65rem; color: #64748b; letter-spacing: 0.05em; }}
            .hist-container {{ margin-top: 25px; background: #f8fafc; padding: 20px; border: 1px solid #e2e8f0; border-radius: 8px; }}
            .summary-overview {{ margin-top: 30px; }}
            .overview-table {{ width: 100%; border-collapse: collapse; font-size: 0.85rem; background: white; border-radius: 8px; overflow: hidden; }}
            .overview-table th {{ background: #f1f5f9; text-align: left; padding: 12px; color: #475569; }}
            .overview-table td {{ padding: 10px 12px; border-bottom: 1px solid #f1f5f9; }}
            .card {{ background: var(--card); border: 1px solid #e2e8f0; border-radius: 12px; margin-bottom: 60px; overflow: hidden; box-shadow: 0 4px 6px -1px rgba(0,0,0,0.05); }}
            .card-header {{ background: #fdfdfd; padding: 15px 25px; border-bottom: 1px solid #f1f5f9; display: flex; justify-content: space-between; align-items: center; }}
            .badge {{ background: var(--accent); color: white; padding: 4px 12px; border-radius: 999px; font-weight: 600; font-size: 0.75rem; }}
            .size-label {{ color: #64748b; font-size: 0.85rem; font-weight: 500; }}
            .card-body {{ display: grid; grid-template-columns: 1fr 1.2fr; gap: 30px; padding: 30px; }}
            .module-title {{ margin: 0 0 15px 0; color: var(--primary); font-size: 1.8rem; border-left: 4px solid var(--accent); padding-left: 15px; }}
            .summary {{ font-size: 1.05rem; font-weight: 600; color: #334155; }}
            .desc {{ font-size: 0.95rem; line-height: 1.7; color: #475569; }}
            .metadata-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 15px; margin: 20px 0; background: #f8fafc; padding: 15px; border-radius: 8px; font-size: 0.85rem; }}
            .pathway-list {{ columns: 1; padding-left: 20px; margin-top: 10px; color: #475569; font-size: 0.9rem; }}
            .string-img {{ width: 100%; border-radius: 8px; border: 1px solid #e2e8f0; }}
            .no-viz {{ background: #f1f5f9; border: 2px dashed #cbd5e1; height: 350px; display: flex; align-items: center; justify-content: center; color: #94a3b8; font-style: italic; border-radius: 8px; }}
            .gene-footer {{ background: #f8fafc; padding: 15px 25px; font-family: monospace; font-size: 0.75rem; color: #64748b; border-top: 1px solid #f1f5f9; }}
            h1 {{ font-size: 2.2rem; margin-top: 0; color: var(--primary); margin-bottom: 10px; }}
        </style>
    </head>
    <body>
    <nav>
        <h3>Module Navigator</h3>
        {nav_links}
    </nav>
    <main>
        <div class="dashboard-header">
            <h1>Unsupervised Module Annotation</h1>
            <div class="subtitle">Generated by Gemini v1.5 Pro | Please Verify Biological Insights</div>
            <div class="stats">
                <div class="stat-box">
                    <div class="stat-val">{len(df)}</div>
                    <div class="stat-label">Total Clusters</div>
                </div>
                <div class="stat-box">
                    <div class="stat-val">{unique_genes_count}</div>
                    <div class="stat-label">Unique Genes</div>
                </div>
                <div class="stat-box">
                    <div class="stat-val">{avg_size}</div>
                    <div class="stat-label">Avg Size</div>
                </div>
            </div>
            <div class="hist-container">
                <div class="stat-label" style="margin-bottom:15px; font-weight:bold;">Module Size Distribution</div>
                {hist_svg}
            </div>
            <div class="summary-overview">
                <h3>Module Summary Overview</h3>
                <table class="overview-table">
                    <thead><tr><th>ID</th><th>Function/Title</th><th>Size</th><th>Primary Theme</th></tr></thead>
                    <tbody>{summary_rows}</tbody>
                </table>
            </div>
        </div>
        {module_cards}
    </main>
    </body>
    </html>
    """
    
    with open(output_html, 'w') as f:
        f.write(full_html)
    print(f"Report generated: {output_html}")

if __name__ == "__main__":
    # Example usage:
    # metadata = { ... } # The dictionary used in the session
    # generate_module_report('us_module_gene_lists.csv', 'new_report.html', metadata)
    pass
