from stardist.models import StarDist2D
import numpy as np

def load_model(model_name="2D_versatile_fluo"):
    return StarDist2D.from_pretrained(model_name)

def predict_crop(model, crop_array):
    labels, details = model.predict_instances(crop_array)
    probs  = details["prob"].tolist() if "prob" in details else []
    coords = details["coord"].tolist() if "coord" in details else []
    return {"probs": probs, "coords": coords, "n_objects": int(labels.max())}

def predict_full(model, img_flat, nrows, ncols):
    img = np.array(img_flat, dtype=np.float32).reshape(nrows, ncols, order="F")
    img_max = img.max()
    if img_max > 0:
        img = img / img_max
    # Lower prob_thresh to 0.1 to catch blebbing/irregular nuclei
    # nms_thresh lowered to 0.2 to avoid merging touching nuclei
    labels, details = model.predict_instances(
        img,
        n_tiles=(4, 4),
        prob_thresh=0.10,
        nms_thresh=0.20
    )
    labels_flat = labels.flatten(order="F").tolist()
    return {"labels_flat": labels_flat, "probs": details["prob"].tolist() if "prob" in details else [], "n_objects": int(labels.max())}

def get_available_keys(model, crop_array):
    labels, details = model.predict_instances(crop_array)
    return list(details.keys())
