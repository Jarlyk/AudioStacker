// anchorp A utility object that automatically scripts the presentation_rect of objects connected to its outlet to follow anchoring relative to the patcher window size

#include "ext.h"
#include "ext_obex.h"
#include "ext_common.h"
#include "jpatcher_api.h"    // jpatcher_api.h must come before z_dsp.h
#include "jgraphics.h"
#include "z_dsp.h"
#include <math.h>
#include "ext_boxstyle.h"

////////////////////////// object struct
typedef struct {
	t_object a_obj;

	int a_anchor_left;
	int a_anchor_right;
	int a_anchor_top;
	int a_anchor_bottom;
	int a_rzoom_left;
	int a_rzoom_right;
	int a_rzoom_top;
	int a_rzoom_bottom;
	
	t_object* a_patcher;
	t_object* a_patcherview;
	t_object* a_box;
	
	void* a_outlet;

	double a_width;
	double a_height;
	double a_zoom;
} t_anchorp;

///////////////////////// function prototypes

void* anchorp_new(t_symbol* s, long argc, t_atom* argv);
void  anchorp_free(t_anchorp* x);
void  anchorp_assist(t_anchorp* x, void* b, long m, long a, char* s);
void  anchorp_notify(t_anchorp* x, t_symbol* s, t_symbol* msg, void* sender, void* data);
void  anchorp_attach(t_anchorp *x);

//////////////////////// global class pointer variable
void* anchorp_class;

void ext_main(void* r) {
	//common_symbols_init();
	t_class* c;
	c = class_new("anchorp", (method)anchorp_new, (method)anchorp_free, sizeof(t_anchorp), (method)NULL, A_GIMME, 0L);

	class_addmethod(c, (method)anchorp_assist, "assist", A_CANT, 0);
	class_addmethod(c, (method)anchorp_notify, "notify", A_CANT, 0);

	CLASS_ATTR_INT32(c, "anchor_left", 0, t_anchorp, a_anchor_left);
	CLASS_ATTR_INT32(c, "anchor_right", 0, t_anchorp, a_anchor_right);
	CLASS_ATTR_INT32(c, "anchor_top", 0, t_anchorp, a_anchor_top);
	CLASS_ATTR_INT32(c, "anchor_bottom", 0, t_anchorp, a_anchor_bottom);
	CLASS_ATTR_INT32(c, "rzoom_left", 0, t_anchorp, a_rzoom_left);
	CLASS_ATTR_INT32(c, "rzoom_right", 0, t_anchorp, a_rzoom_right);
	CLASS_ATTR_INT32(c, "rzoom_top", 0, t_anchorp, a_rzoom_top);
	CLASS_ATTR_INT32(c, "rzoom_bottom", 0, t_anchorp, a_rzoom_bottom);
	CLASS_ATTR_STYLE_LABEL(c, "anchor_left", 0, "onoff", "Anchor Left");
	CLASS_ATTR_STYLE_LABEL(c, "anchor_right", 0, "onoff", "Anchor Right");
	CLASS_ATTR_STYLE_LABEL(c, "anchor_top", 0, "onoff", "Anchor Top");
	CLASS_ATTR_STYLE_LABEL(c, "anchor_bottom", 0, "onoff", "Anchor Bottom");
	CLASS_ATTR_STYLE_LABEL(c, "rzoom_left", 0, "onoff", "Resist Zoom Left");
	CLASS_ATTR_STYLE_LABEL(c, "rzoom_right", 0, "onoff", "Resist Zoom Right");
	CLASS_ATTR_STYLE_LABEL(c, "rzoom_top", 0, "onoff", "Resist Zoom Top");
	CLASS_ATTR_STYLE_LABEL(c, "rzoom_bottom", 0, "onoff", "Resist Zoom Bottom");
	CLASS_ATTR_SAVE(c, "anchor_left", 0);
	CLASS_ATTR_SAVE(c, "anchor_right", 0);
	CLASS_ATTR_SAVE(c, "anchor_top", 0);
	CLASS_ATTR_SAVE(c, "anchor_bottom", 0);
	CLASS_ATTR_SAVE(c, "rzoom_left", 0);
	CLASS_ATTR_SAVE(c, "rzoom_right", 0);
	CLASS_ATTR_SAVE(c, "rzoom_top", 0);
	CLASS_ATTR_SAVE(c, "rzoom_bottom", 0);
	
	anchorp_class = c;
	class_register(CLASS_BOX, c);
}

void anchorp_assist(t_anchorp* x, void* b, long m, long a, char* s) {
	if (m == ASSIST_INLET) {
	}
	else {
		if (a == 0) {
			strncpy_zero(s, "Anchored Objects", 512);
		}
	}
}

void anchorp_notify(t_anchorp *x, t_symbol *s, t_symbol *msg, void *sender, void *data)
{
	if (msg == gensym("attr_modified") && sender == x->a_patcherview) { // sent when an attribute of the sender changes
		t_symbol *attrname;

		// 'data' arg is the modified attribute object
		// get its name:
		attrname = (t_symbol *)object_method(data, gensym("getname"));

		if (attrname == gensym("rect")) {
			t_atom *av = NULL;
			long ac = 0;

			object_attr_getvalueof(sender, attrname, &ac, &av);
			if (ac && av) {
				double new_width = atom_getfloat(av+2);
				double new_height = atom_getfloat(av+3);
				if (new_width != x->a_width || new_height != x->a_height) {

					//Enumerate patch cords, looking for those originating from our box
					t_object* line = jpatcher_get_firstline(x->a_patcher);
					while (line != NULL) {
						if (jpatchline_get_box1(line) == x->a_box) {
							t_object* box2 = jpatchline_get_box2(line);
							t_rect pr;
							if (jbox_get_presentation_rect(box2, &pr) == MAX_ERR_NONE) {
								double dx = new_width - x->a_width;
								double dy = new_height - x->a_height;

								if (!x->a_anchor_left){
									pr.x += dx;
								}
								if (!x->a_anchor_top) {
									pr.y += dy;
								}
								if (x->a_anchor_right && x->a_anchor_left) {
									pr.width += dx;
								}
								if (x->a_anchor_bottom && x->a_anchor_top) {
									pr.height += dy;
								}
								jbox_set_presentation_rect(box2, &pr);									
							}							
						}
						line = jpatchline_get_nextline(line);
					}

					x->a_width = new_width;
					x->a_height = new_height;
				}
				freebytes(av, sizeof(t_atom) * ac);
			}
		}
		else if (attrname == gensym("zoomfactor")) {
			double new_zoom = patcherview_get_zoomfactor(x->a_patcherview);		
			if (new_zoom != x->a_zoom) {
				//We're going to need the patcher size to perform the zoom rescaling
				t_rect rectp;
				patcherview_get_rect(x->a_patcherview, &rectp);
				
				//Enumerate patch cords, looking for those originating from our box
				t_object* line = jpatcher_get_firstline(x->a_patcher);
				while (line != NULL) {
					if (jpatchline_get_box1(line) == x->a_box) {
						t_object* box2 = jpatchline_get_box2(line);
						t_rect pr;
						if (jbox_get_presentation_rect(box2, &pr) == MAX_ERR_NONE) {
							//Our goal when resisting a zoom change is to maintain the margin size the same in displayed pixels
							//This means we want to get larger as zoom gets smaller and vice versa
							if (x->a_rzoom_left){
								pr.x *= (x->a_zoom/new_zoom);
							}
							if (x->a_rzoom_top) {
								pr.y *= (x->a_zoom/new_zoom);
							}
							if (x->a_rzoom_right) {
								pr.width *= (x->a_zoom/new_zoom);
							}
							if (x->a_rzoom_bottom) {
								pr.height *= (x->a_zoom/new_zoom);
							}
							jbox_set_presentation_rect(box2, &pr);									
						}							
					}
					line = jpatchline_get_nextline(line);
				}
				x->a_zoom = new_zoom;

				if (x->a_rzoom_left || x->a_rzoom_right || x->a_rzoom_top || x->a_rzoom_bottom) {
					//Try to disable scroll bars whenever zoom changes
					//This resolves an issue where the presentation rects update too late and scroll bars automatically appear					
					object_attr_setchar(x->a_patcher, gensym("enablehscroll"), 0);
					object_attr_setchar(x->a_patcher, gensym("enablevscroll"), 0);					
				}
			}
		}
	}
}

void anchorp_free(t_anchorp* x) {
	t_object *jp = NULL;
	t_object *pv;

	// detach from any objects that you have attached to.
	object_obex_lookup(x, gensym("#P"), &jp);
	if (jp) {
		pv = jpatcher_get_firstview(jp);
		object_detach_byptr(x, pv);
	}
}

void* anchorp_new(t_symbol* s, long argc, t_atom* argv) {
	t_anchorp* x = NULL;

	if ((x = (t_anchorp*)object_alloc(anchorp_class))) {
		x->a_outlet = outlet_new(x, NULL);

		//Default anchoring doesn't manipulate the views
		x->a_anchor_left = 1;
		x->a_anchor_right = 0;
		x->a_anchor_top = 1;
		x->a_anchor_bottom = 0;
		x->a_rzoom_left = 0;
		x->a_rzoom_right = 0;
		x->a_rzoom_top = 0;
		x->a_rzoom_bottom = 0;
		
		// Get a pointer to our patcher and box
		object_obex_lookup(x, gensym("#P"), &x->a_patcher);

		// The patcherview is not available when the object is created as a patcher is being read from disk,
		// so we have to defer to wait for it before getting it and attaching for notifications.
		defer_low(x, (method)anchorp_attach, NULL, 0, NULL);
	}
	
	return x;
}

void anchorp_attach(t_anchorp *x)
{
	t_atom *av = NULL;
	long ac = 0;

	object_obex_lookup(x, gensym("#B"), &x->a_box);
	
	x->a_patcherview = object_attr_getobj(x->a_patcher, gensym("firstview"));
	object_attach_byptr_register(x, x->a_patcherview, CLASS_NOBOX);

	// get the bounds of the first patcherview and cache them
	object_attr_getvalueof(x->a_patcherview, gensym("rect"), &ac, &av);
	if (ac && av) {
		x->a_width = atom_getfloat(av+2);
		x->a_height = atom_getfloat(av+3);
		freebytes(av, sizeof(t_atom) * ac); // or sysmem_freeptr()
	}
	x->a_zoom = patcherview_get_zoomfactor(x->a_patcherview);
}