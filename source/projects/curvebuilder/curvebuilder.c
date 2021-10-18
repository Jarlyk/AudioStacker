// stackscope~ A triggered time-domain scope that shows relative contributions of two signals to a sum, with support for stereo channels

#include "ext.h"
#include "ext_obex.h"
#include "ext_common.h"
#include "jpatcher_api.h"    // jpatcher_api.h must come before z_dsp.h
#include "jgraphics.h"
#include "z_dsp.h"
#include <math.h>
#include "ext_boxstyle.h"
#include "utils.h"

////////////////////////// object struct
typedef struct {
	t_pxjbox a_obj;

	//Flag that's set true if an internal error breaks the scope to the point that it can no longer function
	int a_faulted;

	//Whether currently active (when this is disabled, will ignore signals and only draw the background)
	int a_active;
	
	//Background color
	t_jrgba a_bgcolor;
	t_jrgba a_sig1color;
	t_jrgba a_sigsumcolor;
	t_jrgba a_aboveceilcolor;

	//Display range of the chart, in dB
	double a_range;
	double a_range_s;

	//Ceiling indicator for chart, in dB; samples above this threshold will be highlighted in red
	double a_ceiling;
	double a_ceiling_s;

	//Which channel to display (mono, left, right, max)
	t_symbol* a_channel;
	
	//Raw signal buffers with current snapshot
	//Note that these can be updated live and thus have a mixture of prior and new snapshot
	double* a_sig1_left;
	double* a_sig1_right;
	double* a_sig2_left;
	double* a_sig2_right;

	//Current sampling index into signal buffers
	int a_sig_index;

	//Request for new sample count
	int a_pending_samples;
		
	//Buffer for rendering of chart
	t_uint8* a_img_chart;

	//Downsampled buffer for AA
	t_uint8* a_img_chart_ds;

	//Current image size of full buffer (twice as large as the downsampled buffer that's actually displayed)
	int   a_sizex;
	int   a_sizey;

	//Number of samples to capture per snapshot
	int   a_samples;

	//Whether currently mid-sampling
	bool  a_is_sampling;

	//Bang outlet when a snapshot was fully captured
	void *a_outlet1_snapped;

	//Float outlets with peak of left/right sum signals, in dB
	void *a_outlet2_peak_left;
	void *a_outlet3_peak_right;
} t_stackscope;

///////////////////////// function prototypes

void stackscope_dsp64(t_stackscope* x, t_object* dsp64, short* count, double samplerate, long maxvectorsize, long flags);
void stackscope_perform64(t_stackscope* x, t_object* dsp64, double** ins, long numins, double** outs, long numouts, long sampleframes,
	long flags, void* userparam);

void stackscope_paint(t_stackscope* x, t_object* view);

void* stackscope_new(t_symbol* s, long argc, t_atom* argv);
void  stackscope_free(t_stackscope* x);
void  stackscope_assist(t_stackscope* x, void* b, long m, long a, char* s);

// Helper functions
void update_img_chart(t_stackscope* scope, int sizex, int sizey);
void draw_column_with_ceil(t_stackscope* scope, t_uint8* buf, int x, int stride, double v0, double v1, t_jrgba c);
void draw_column(t_stackscope* scope, t_uint8* buf, int x, int stride, double v0, double v1, t_jrgba c);
int gety(t_stackscope* scope, double v);

//////////////////////// global class pointer variable
void* stackscope_class;

//Other static globals
t_symbol* sym_mono;
t_symbol* sym_left;
t_symbol* sym_right;
t_symbol* sym_max;

void ext_main(void* r) {
	sym_mono = gensym("mono");
	sym_left = gensym("left");
	sym_right = gensym("right");
	sym_max = gensym("max");
	
	//common_symbols_init();
	t_class* c;
	c = class_new("stackscope~", (method)stackscope_new, (method)stackscope_free, sizeof(t_stackscope), (method)NULL, A_GIMME, 0L);

	c->c_flags |= CLASS_FLAG_NEWDICTIONARY; /* | CLASS_FLAG_VISUALIZER; */
	jbox_initclass(c, 0);
	class_dspinitjbox(c);

	class_addmethod(c, (method)stackscope_dsp64, "dsp64", A_CANT, 0);
	class_addmethod(c, (method)stackscope_paint, "paint", A_CANT, 0);
	class_addmethod(c, (method)stackscope_assist, "assist", A_CANT, 0);
	class_addmethod(c, (method)jbox_notify, "notify", A_CANT, 0);

	CLASS_ATTR_INT32(c, "active", 0, t_stackscope, a_active);
	CLASS_ATTR_RGBA(c, "bgcolor", 0, t_stackscope, a_bgcolor);
	CLASS_ATTR_RGBA(c, "sig1color", 0, t_stackscope, a_sig1color);
	CLASS_ATTR_RGBA(c, "sigsumcolor", 0, t_stackscope, a_sigsumcolor);
	CLASS_ATTR_RGBA(c, "aboveceilcolor", 0, t_stackscope, a_aboveceilcolor);
	CLASS_ATTR_INT32(c, "samples", 0, t_stackscope, a_pending_samples);
	CLASS_ATTR_FILTER_MIN(c, "samples", 0);
	CLASS_ATTR_DOUBLE(c, "range", 0, t_stackscope, a_range);
	CLASS_ATTR_DOUBLE(c, "ceiling", 0, t_stackscope, a_ceiling);

	CLASS_ATTR_SYM(c, "channel", 0, t_stackscope, a_channel);
	CLASS_ATTR_ENUM(c, "channel", 0, "mono left right max");
	CLASS_ATTR_DEFAULT(c, "channel", 0, "mono");

	CLASS_ATTR_DEFAULT(c, "patching_rect", 0, "0. 0. 120. 80.");

	stackscope_class = c;
	class_register(CLASS_BOX, c);

}

void stackscope_dsp64(t_stackscope* x, t_object* dsp64, short* count, double samplerate, long maxvectorsize, long flags) {
	dsp_add64(dsp64, (t_object*)x, (t_perfroutine64)stackscope_perform64, 0, NULL);
}


void stackscope_perform64(t_stackscope* x, t_object* dsp64, double** ins, long numins, double** outs, long numouts, long sampleframes,
	long flags, void* userparam) {
	
	if (x->a_faulted) return;
	
	if (numins != 5) {
		object_error((t_object*)x, "Invalid number of signal inputs into stackscope_perform64");		
		x->a_faulted = true;
		return;
	}

	double* sig1_l = ins[0];
	double* sig1_r = ins[1];
	double* sig2_l = ins[2];
	double* sig2_r = ins[3];
	double* trigger = ins[4];
	int search_idx = 0;
	int samp_idx = 0;

	while (search_idx < sampleframes) {
		if (x->a_is_sampling) {
			int new_index = x->a_sig_index + (sampleframes - samp_idx);
			if (new_index >= x->a_samples) {
				int n = x->a_samples - x->a_sig_index;
				size_t size = n*sizeof(double);
				memcpy(&x->a_sig1_left[x->a_sig_index], &sig1_l[samp_idx], size);
				memcpy(&x->a_sig1_right[x->a_sig_index], &sig1_r[samp_idx], size);
				memcpy(&x->a_sig2_left[x->a_sig_index], &sig2_l[samp_idx], size);
				memcpy(&x->a_sig2_right[x->a_sig_index], &sig2_r[samp_idx], size);
				x->a_sig_index = x->a_samples;
				x->a_is_sampling = false;
				search_idx = samp_idx + n;

				//Report peaks and snap completion
				double peak_left_s = 0;
				double peak_right_s = 0;
				for (int i=0; i < x->a_samples; i++) {
					peak_left_s = max(peak_left_s, fabs(x->a_sig1_left[i] + x->a_sig2_left[i]));
					peak_right_s = max(peak_right_s, fabs(x->a_sig1_right[i] + x->a_sig2_right[i]));
				}
				outlet_float(x->a_outlet2_peak_left, atodb(peak_left_s));
				outlet_float(x->a_outlet3_peak_right, atodb(peak_right_s));
				outlet_bang(x->a_outlet1_snapped);
			}
			else {
				size_t size = (sampleframes - samp_idx)*sizeof(double);
				memcpy(&x->a_sig1_left[x->a_sig_index], &sig1_l[samp_idx], size);
				memcpy(&x->a_sig1_right[x->a_sig_index], &sig1_r[samp_idx], size);
				memcpy(&x->a_sig2_left[x->a_sig_index], &sig2_l[samp_idx], size);
				memcpy(&x->a_sig2_right[x->a_sig_index], &sig2_r[samp_idx], size);
				x->a_sig_index += (sampleframes - samp_idx);
				search_idx = sampleframes;
			}

			jbox_redraw((t_jbox*)x);
		}

		if (!x->a_is_sampling) {
			while (search_idx < sampleframes) {
				if (trigger[search_idx] > 0) {
					if (x->a_pending_samples > 0 && (x->a_samples != x->a_pending_samples || x->a_sig1_left == NULL))
					{
						x->a_samples = x->a_pending_samples;
						if (x->a_sig1_left != NULL) sysmem_freeptr(x->a_sig1_left);
						if (x->a_sig1_right != NULL) sysmem_freeptr(x->a_sig1_right);
						if (x->a_sig2_left != NULL) sysmem_freeptr(x->a_sig2_left);
						if (x->a_sig2_right != NULL) sysmem_freeptr(x->a_sig2_right);
						x->a_sig1_left = (double*)sysmem_newptr(x->a_samples*sizeof(double));
						x->a_sig1_right = (double*)sysmem_newptr(x->a_samples*sizeof(double));
						x->a_sig2_left = (double*)sysmem_newptr(x->a_samples*sizeof(double));
						x->a_sig2_right = (double*)sysmem_newptr(x->a_samples*sizeof(double));
					}
					
					x->a_sig_index = 0;
					samp_idx = search_idx;
					x->a_is_sampling = true;
					break;
				}
				search_idx++;
			}
		}		
	}
}

void stackscope_paint(t_stackscope* x, t_object* view) {
	if (!x->a_active) return;
	
	t_jgraphics* g;
	t_rect       rect;

	//Set up graphics interface
	g = (t_jgraphics*)patcherview_get_jgraphics(view);	
	jbox_get_rect_for_view((t_object*)x, view, &rect);
    double zoom = patcherview_get_zoomfactor(view);	
	rect.x = rect.y = 0;

	//Draw background rectangle
	if (x->a_bgcolor.alpha > 0) {
		jgraphics_set_source_jrgba(g, &x->a_bgcolor);
		jgraphics_rectangle_fill_fast(g, 0, 0, rect.width, rect.height);		
	}

	//Update chart image based on current size
	update_img_chart(x, (int)(2*zoom*rect.width), (int)(2*zoom*rect.height));
	
	if (x->a_sizex <= 0 || x->a_sizey <= 0 || x->a_img_chart_ds == NULL) return;
	
	//Draw chart image
	t_jsurface* surf = jgraphics_image_surface_create_for_data(x->a_img_chart_ds, 
		JGRAPHICS_FORMAT_ARGB32,
		x->a_sizex/2,
		x->a_sizey/2,
		2*x->a_sizex,
		NULL,
		NULL);
	jgraphics_image_surface_draw_fast(g, surf);
	jgraphics_surface_destroy(surf);

	//TODO: Draw annotation layer over the chart image
	//This can include things like mouseover time/dB display
}

void draw_full_column(int sizey, int stride, t_uint8* bufIter, t_uint8 b, t_uint8 g, t_uint8 r, t_uint8 a) {
	for (int y=0; y < sizey; y++) {
		bufIter[0] = b;
		bufIter[1] = g;
		bufIter[2] = r;
		bufIter[3] = a;
		bufIter += stride;
	}
}

typedef double (*f_get_channel)(double, double);

double get_mono(double left, double right) {
	return 0.5*(left + right);
}

double get_left(double left, double right) {
	return left;
}

double get_right(double left, double right ) {
	return right;
}

double get_max_lr(double left, double right) {
	return max(left, right);
}

f_get_channel get_channel_map(t_stackscope* scope) {
	if (scope->a_channel == sym_left) return &get_left;
	if (scope->a_channel == sym_right) return &get_right;
	if (scope->a_channel == sym_max) return &get_max_lr;
	return &get_mono;
}

void update_img_chart(t_stackscope* scope, int sizex, int sizey) {
	if (sizex <= 0 || sizey <= 0) return;
	
	if (sizex != scope->a_sizex || sizey != scope->a_sizey || scope->a_img_chart == NULL) {
		if (scope->a_img_chart != NULL) {
			sysmem_freeptr(scope->a_img_chart);
		}
		if (scope->a_img_chart_ds != NULL) {
			sysmem_freeptr(scope->a_img_chart_ds);
		}
		scope->a_sizex = sizex;
		scope->a_sizey = sizey;
		scope->a_img_chart = (t_uint8*)sysmem_newptrclear(4*sizex*sizey);
		scope->a_img_chart_ds = (t_uint8*)sysmem_newptrclear(sizex*sizey);
	}

	scope->a_range_s = dbtoa(scope->a_range);
	scope->a_ceiling_s = dbtoa(scope->a_ceiling);

	//TODO: Support non-full-redraw

	//Create some local helper variables to keep things concise
	t_uint8* buf = scope->a_img_chart;
	int nsig = scope->a_samples;
	int sic = scope->a_sig_index;
	double* s1l = scope->a_sig1_left;
	double* s1r = scope->a_sig1_right;
	double* s2l = scope->a_sig2_left;
	double* s2r = scope->a_sig2_right;

	//Create channel mapper function
	f_get_channel cmap = get_channel_map(scope);

	//Compute a few things we'll need
	int lastx = sic*(sizex-1)/(nsig-1);
	if (lastx >= sizex) lastx = sizex-1;
	int stride = 4*sizex;

	//Create image one column at a time
	double ext1_prior = 0;
	double ext2_prior = 0;
	double extsum_prior = 0;
	for (int x=0; x <= lastx; x++) {
		
		//Start by marking all pixels in the column as transparent
		t_uint8* bufIter = buf + 4*x;
		draw_full_column(sizey, stride, bufIter, 0, 0, 0, 0);
		
		//Get first and last sample index for this column's bin
		int si0 = x*(nsig-1)/(sizex-1);
		int si1 = (x+1)*(nsig-1)/(sizex-1) - 1;
		if (si1 < si0) si1 = si0 + 1;
		if (si1 >= nsig) si1 = nsig - 1;		

		//Get min/max values for this bin (TODO: support other downsampling methods)
		double min1 = DBL_MAX;
		double max1 = DBL_MIN;
		double min2 = DBL_MAX;
		double max2 = DBL_MIN;
		double minsum = DBL_MAX;
		double maxsum = DBL_MIN;
		for (int i=si0; i<=si1; i++) {
			double s1 = cmap(s1l[i], s1r[i]);
			double s2 = cmap(s2l[i], s2r[i]);
			min1 = min(s1, min1);
			max1 = max(s1, max1);
			min2 = min(s2, min2);
			max2 = max(s2, max2);
			minsum = min(s1 + s2, minsum);
			maxsum = max(s1 + s2, maxsum);
		}

		//Compute extreme edges
		double ext1 = max1 > -1*min1 ? max1 : min1;
		double ext2 = max2 > -1*min2 ? max2 : min2;		
		double extsum = maxsum > -1*minsum ? maxsum : minsum;

		//We only start drawing after the first bin is computed
		//This allows us to connect large y displacements
		if (x > 0) {
			double ext1_mid = 0.5*(ext1_prior + ext1);
			draw_column(scope, buf, x-1, stride, ext1_prior, ext1_mid, scope->a_sig1color);
			draw_column(scope, buf, x, stride, ext1_mid, ext1, scope->a_sig1color);
			
			double extsum_mid = 0.5*(extsum_prior + extsum);
			draw_column_with_ceil(scope, buf, x-1, stride, extsum_prior, extsum_mid, scope->a_sigsumcolor);
			draw_column_with_ceil(scope, buf, x, stride, extsum_mid, extsum, scope->a_sigsumcolor);
		}

		ext1_prior = ext1;
		ext2_prior = ext2;
		extsum_prior = extsum;
	}

	//Draw leading edge indicator if there's room
	if (lastx + 1 < sizex)
		draw_full_column(sizey, stride, buf + 4*(lastx + 1), 0, 255, 0, 255);
	if (lastx + 2 < sizex)
		draw_full_column(sizey, stride, buf + 4*(lastx + 2), 0, 128, 0, 255);
	
	//Now downsample to anti-alias
	t_uint8* bufIter = buf;
	t_uint8* dsIter = scope->a_img_chart_ds;
	for (int y=0; y < sizey/2; y++) {
		for (int x = 0; x < sizex/2; x++) {
			for (int k=0; k < 4; k++) {
				*dsIter = (t_uint8)(((int)bufIter[0] + bufIter[4] + bufIter[stride] + bufIter[stride+4])/4);
				dsIter++;
				bufIter++;				
			}
			
			bufIter += 4;
		}

		bufIter += 4*sizex;
	}
}

void draw_pixel(t_uint8* bufIter, t_jrgba c) {
	//Compute color from interpolating along gradient
	bufIter[0] = (t_uint8)floor(255*c.blue);
	bufIter[1] = (t_uint8)floor(255*c.green);
	bufIter[2] = (t_uint8)floor(255*c.red);
	bufIter[3] = (t_uint8)floor(255*c.alpha);
}

void draw_column_with_ceil(t_stackscope* scope, t_uint8* buf, int x, int stride, double v0, double v1, t_jrgba c) {
	int v0greater = fabs(v0) > scope->a_ceiling_s;
	int v1greater = fabs(v1) > scope->a_ceiling_s;
	if (v0greater && v1greater) {
		draw_column(scope, buf, x, stride, v0, v1, scope->a_aboveceilcolor);
	}
	else if (v0greater) {
		draw_column(scope, buf, x, stride, v1, scope->a_ceiling_s, c);
		draw_column(scope, buf, x, stride, v0, scope->a_ceiling_s, scope->a_aboveceilcolor);		
	}
	else if (v1greater) {
		draw_column(scope, buf, x, stride, v0, scope->a_ceiling_s, c);
		draw_column(scope, buf, x, stride, v1, scope->a_ceiling_s, scope->a_aboveceilcolor);
	}
	else {
		draw_column(scope, buf, x, stride, v0, v1, c);
	}
}

void draw_column(t_stackscope* scope, t_uint8* buf, int x, int stride, double v0, double v1, t_jrgba c) {
	int y0 = gety(scope, v0);
	int y1 = gety(scope, v1);

	if (y0 <= y1) {
		t_uint8* bufIter = buf + y0*stride + 4*x;
		for (int y=y0; y <= y1; y++) {
			draw_pixel(bufIter, c);
			bufIter += stride;
		}
	}
	else {
		t_uint8* bufIter = buf + y1*stride + 4*x;
		for (int y=y1; y <= y0; y++) {
			draw_pixel(bufIter, c);
			bufIter += stride;
		}		
	}
}

int gety(t_stackscope* scope, double v) {
	int sizey = scope->a_sizey;
	int y = (int)round((0.5 - 0.5*v/scope->a_range_s)*sizey);
	if (y < 0) y = 0;
	if (y >= sizey) y = sizey - 1;
	return y;
}

void stackscope_assist(t_stackscope* x, void* b, long m, long a, char* s) {
	if (m == ASSIST_INLET) {
		switch (a) {
			case 0:
				strncpy_zero(s, "Signal 1 (L)", 512);
				break;
			case 1:
				strncpy_zero(s, "Signal 1 (R)", 512);
				break;
			case 2:
				strncpy_zero(s, "Signal 2 (L)", 512);
				break;
			case 3:
				strncpy_zero(s, "Signal 2 (R)", 512);
				break;
			case 4:
				strncpy_zero(s, "Trigger", 512);
				break;
		}
	}
	else {
	}
}

void stackscope_free(t_stackscope* x) {
	if (x->a_img_chart != NULL) {
		sysmem_freeptr(x->a_img_chart);
		x->a_img_chart = NULL;
	}
	if (x->a_img_chart_ds != NULL) {
		sysmem_freeptr(x->a_img_chart_ds);
		x->a_img_chart_ds = NULL;
	}
	if (x->a_sig1_left != NULL) {
		sysmem_freeptr(x->a_sig1_left);
		x->a_sig1_left = NULL;
	}
	if (x->a_sig1_right != NULL) {
		sysmem_freeptr(x->a_sig1_right);
		x->a_sig1_right = NULL;
	}
	if (x->a_sig2_left != NULL) {
		sysmem_freeptr(x->a_sig2_left);
		x->a_sig2_left = NULL;
	}
	if (x->a_sig2_right != NULL) {
		sysmem_freeptr(x->a_sig2_right);
		x->a_sig2_right = NULL;
	}
	
	dsp_freejbox((t_pxjbox*)x);
	jbox_free((t_jbox*)x);
}

void* stackscope_new(t_symbol* s, long argc, t_atom* argv) {
	t_stackscope* x = NULL;
	t_dictionary* d = NULL;

	if (!(d = object_dictionaryarg(argc, argv)))
		return NULL;
	long i;

	if ((x = (t_stackscope*)object_alloc(stackscope_class))) {
		long boxflags = 0 | JBOX_DRAWFIRSTIN | JBOX_NODRAWBOX | JBOX_DRAWINLAST | JBOX_GROWBOTH | JBOX_DRAWBACKGROUND;

		jbox_new((t_jbox*)x, boxflags, argc, argv);

		//Default values
		x->a_sig_index = 0;
		jrgba_set(&x->a_bgcolor, 0, 0, 0, 1);
		jrgba_set(&x->a_sig1color, 0, 0.3, 1, 1);
		jrgba_set(&x->a_sigsumcolor, 0 ,1, 0, 1);
		jrgba_set(&x->a_aboveceilcolor, 1, 0, 0, 1);
		x->a_active = 1;
		x->a_range = 1;
		x->a_ceiling = 0;
		x->a_img_chart = NULL;
		x->a_is_sampling = false;
		x->a_pending_samples = 4800;
		x->a_samples = 0;
		x->a_sig1_left = NULL;
		x->a_sig1_right = NULL;
		x->a_sig2_left = NULL;
		x->a_sig2_right = NULL;
		x->a_sizex = 0;
		x->a_sizey = 0;
		x->a_faulted = false;
		
		x->a_obj.z_box.b_firstin = (void*)x;
		dsp_setupjbox((t_pxjbox*)x, 5);

		x->a_outlet1_snapped = bangout((t_object *)x);
		x->a_outlet2_peak_left = floatout((t_object *)x);
		x->a_outlet3_peak_right = floatout((t_object *)x);

		attr_dictionary_process(x, d);
		jbox_ready((t_jbox*)x);
	}
	return x;
}
