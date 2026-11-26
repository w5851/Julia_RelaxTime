"""Read the five preselected Fig. 5-3 features from native PDF paths.

No interpolation or smoothing is applied. Coordinates describe the published
polyline, not the author's original unrounded numerical data. The overlay is
temporary review material; the JSON retains calibration and adjacent vertices.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pymupdf
from PIL import Image, ImageDraw


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("pdf", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("overlay", type=Path)
    args = parser.parse_args()
    doc = pymupdf.open(args.pdf)
    page = doc[82]
    drawings = page.get_drawings()
    candidates = [d for d in drawings if d["color"] == (0., 0., 0.)
                  and d["rect"].y0 > 410 and len(d["items"]) > 300]
    if len(candidates) != 1:
        raise ValueError("Expected one native Fig. 5-3 polyline")
    curve = candidates[0]
    if any(it[0] != "l" for it in curve["items"]):
        raise ValueError("This extractor requires the observed straight polyline")
    vertices = [curve["items"][0][1]] + [it[2] for it in curve["items"]]
    x_ticks = []
    y_ticks = []
    for d in drawings:
        if d["color"] is None:
            continue
        for it in d["items"]:
            if it[0] != "l":
                continue
            a, b = it[1:3]
            if (abs(a.x-b.x) < 1e-4 and 3.2 < abs(a.y-b.y) < 3.4
                    and abs(max(a.y,b.y)-603.229) < .01):
                x_ticks.append(a.x)
            if (abs(a.y-b.y) < 1e-4 and 3.2 < abs(a.x-b.x) < 3.4
                    and abs(min(a.x,b.x)-191.867) < .01
                    and 420 < a.y < 605):
                y_ticks.append(a.y)
    x_ticks = sorted(set(x_ticks))
    y_ticks = sorted(set(y_ticks), reverse=True)
    assert len(x_ticks) == 7 and len(y_ticks) == 9
    x_values = np.arange(1000., 1601., 100.)
    y_values = np.arange(9.)
    x_scale = 600/(x_ticks[-1]-x_ticks[0])
    y_scale = 8/(y_ticks[-1]-y_ticks[0])

    def physical(point):
        return [1000+(point.x-x_ticks[0])*x_scale,
                (point.y-y_ticks[0])*y_scale]

    # The initial almost-empty onset bend precedes the five already selected
    # low-density features. Feature order is fixed by the previous readout.
    turns = [i for i in range(1,len(vertices)-1)
             if (vertices[i].x-vertices[i-1].x)*
                (vertices[i+1].x-vertices[i].x) < 0]
    selected = turns[1:6]
    ids = ("first_minimum", "first_maximum", "second_minimum",
           "second_maximum", "main_minimum")
    features = []
    for name, index in zip(ids,selected):
        p = vertices[index]
        features.append({"id":name,"vertex_index":index,
                         "pdf_point":list(p),"muB_MeV":physical(p)[0],
                         "rhoB_over_rho0":physical(p)[1],
                         "adjacent_vertices":[physical(q) for q in vertices[index-1:index+2]],
                         "half_stroke_muB_MeV":curve["width"]*x_scale/2})
    result = {
        "pdf_path":str(args.pdf.resolve()),
        "pdf_sha256":hashlib.sha256(args.pdf.read_bytes()).hexdigest(),
        "pdf_page":83,"printed_page":69,"figure":"5-3",
        "T_MeV":5.,"eB_GeV2":.05,
        "method":"native PDF polyline vertices, calibrated by vector axis ticks",
        "path_sequence":curve["seqno"],"path_vertex_count":len(vertices),
        "line_width_pdf_points":curve["width"],
        "x_ticks_pdf_points":x_ticks,"x_tick_values_MeV":x_values.tolist(),
        "y_ticks_pdf_points":y_ticks,"y_tick_values":y_values.tolist(),
        "x_check_max_residual_MeV":float(max(abs(1000+(x-x_ticks[0])*x_scale-v)
                                                   for x,v in zip(x_ticks,x_values))),
        "y_check_max_residual":float(max(abs((y-y_ticks[0])*y_scale-v)
                                               for y,v in zip(y_ticks,y_values))),
        "features":features,
        "uncertainty":"Stroke and adjacent vertices describe graphical resolution; no original numerical accuracy or statistical interval is inferred."
    }
    # Re-running a coordinate extraction does not discard the separately
    # verified review when the PDF and exact features remain identical.
    if args.output.exists():
        previous=json.loads(args.output.read_text(encoding="utf-8"))
        if (previous.get("pdf_sha256")==result["pdf_sha256"]
                and previous.get("features")==result["features"]
                and "visual_review" in previous):
            result["visual_review"]=previous["visual_review"]
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    clip = pymupdf.Rect(156,409,440,637)
    zoom=4
    pix=page.get_pixmap(matrix=pymupdf.Matrix(zoom,zoom),clip=clip)
    image=Image.frombytes("RGB",[pix.width,pix.height],pix.samples)
    draw=ImageDraw.Draw(image)
    for index,f in enumerate(features,1):
        px=(f["pdf_point"][0]-clip.x0)*zoom
        py=(f["pdf_point"][1]-clip.y0)*zoom
        draw.ellipse((px-6,py-6,px+6,py+6),outline="red",width=2)
        draw.text((px+8,py-8),str(index),fill="red")
    args.overlay.parent.mkdir(parents=True,exist_ok=True)
    image.save(args.overlay)
    print(json.dumps({"features":features,"x_check_max_residual_MeV":result["x_check_max_residual_MeV"]},ensure_ascii=True))


if __name__ == "__main__":
    main()
