from stardist.models import StarDist2D
import numpy as np

def load_model(model_name='2D_versatile_fluo'):
    return StarDist2D.from_pretrained(model_name)

def predict_crop(model, crop_array):
    labels, details = model.predict_instances(crop_array)
    probs  = details['prob'].tolist()  if 'prob'  in details else []
    coords = details['coord'].tolist() if 'coord' in details else []
    n_obj  = int(labels.max())
    return {'probs': probs, 'coords': coords, 'n_objects': n_obj}

def get_available_keys(model, crop_array):
    labels, details = model.predict_instances(crop_array)
    return list(details.keys())
