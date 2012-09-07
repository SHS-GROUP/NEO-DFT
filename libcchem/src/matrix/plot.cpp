#include "plot.hpp"

#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkStructuredGridGeometryFilter.h"

#include "vtkActor.h"
#include "vtkAxesActor.h"
#include "vtkCamera.h"
#include "vtkCaptionActor2D.h"
#include "vtkContourFilter.h"
#include "vtkDataSetMapper.h"
#include "vtkGlyph3D.h"
#include "vtkGlyphSource2D.h"
#include "vtkOutlineFilter.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkPNGWriter.h"
#include "vtkProperty.h"
#include "vtkProperty2D.h"
#include "vtkRectilinearGridGeometryFilter.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkStructuredGridGeometryFilter.h"
#include "vtkScalarBarActor.h"
#include "vtkTextProperty.h"
#include "vtkTubeFilter.h"
#include "vtkWarpScalar.h"
#include "vtkWindowToImageFilter.h"

#include "contrib/matplot/matplot.h"

#include <boost/tuple/tuple.hpp>

namespace matrix {
    namespace plot {

	struct index_vector {
	    index_vector(size_t size, int from = 0) : size_(size), from_(from) {}
	    size_t size() const { return size_; }
	    int operator()(int i) const { return from_ + i; }
	private: size_t size_; int from_;
	};

	class Renderer {
	public:
	    Renderer(size_t xPix, size_t yPix) : started_(false) {
		this->rend = vtkRenderer::New();
		this->grid = vtkStructuredGrid::New();
		this->renWin = vtkRenderWindow::New();
		this->iren = vtkRenderWindowInteractor::New();
		renWin->AddRenderer(this->rend);
		iren->SetRenderWindow(renWin);
		renWin->SetSize(xPix, yPix);

		// Start interactive rendering
		//iren->Initialize();
		//iren->Start();

	    }
	    ~Renderer() {
		iren->Delete();
		renWin->Delete();
		grid->Delete();
		rend->Delete();
	    }
	    void start() {
		//std::cout << "start();\n";
		started_ = true;
		iren->Initialize();
		iren->Start();
	    }
	    bool started() const { return started_; }
	    operator vtkRenderer*() { return rend; }
	    operator vtkStructuredGrid*() { return grid; }
	    void render() {
		renWin->Render();
	    }
	private:
	    bool started_;
	    vtkRenderer *rend;
	    vtkStructuredGrid *grid;
	    vtkRenderWindow *renWin;
	    vtkRenderWindowInteractor *iren;
	};

	template<class V, class M>
	double geometry(const V &x, const V &y, const M &z, vtkStructuredGrid *grid) {
	    const size_t Nx = x.size();
	    const unsigned int Ny = y.size();
	    unsigned int i, j, k;

	    // make sure the input is ok and that this surfaceplot is free
	    assert(Nx == z.size1());
	    assert(Ny == z.size2());

	    // put data, z, into a 2D structured grid
	    grid->SetDimensions(Nx, Ny, 1);

	    double z_low = 10000, z_upp = -10000;

	    vtkPoints *points = vtkPoints::New();
	    points->Allocate(Nx*Ny);
	    for (j = 0, k = 0; j < Ny; j++) {
		for (i = 0; i < Nx; i++, ++k) {
		    //std::cout << x(i)  << " " << y(j) <<  " " << z(i, j) << "\n";
		    points->InsertPoint(k, x(i), y(j), z(i, j));

		    if (z(i, j)< z_low) z_low = z(i, j);
		    if (z(i, j)> z_upp) z_upp = z(i, j);
		}
	    }
	    grid->SetPoints(points);

	    double Lz = z_upp-z_low;
	    // determine x-y range of data
	    double Lxy = std::max(x(Nx-1)-x(0), y(Ny-1)-y(0));
	    double scale = Lxy/Lz;
	    

	    // get scalar field from z-values
	    vtkFloatArray *colors = vtkFloatArray::New();
	    colors->SetNumberOfComponents(1);
	    colors->SetNumberOfTuples(Nx*Ny);
	    for (j = 0, k = 0; j < Ny; j++) {
		for (i = 0; i < Nx; i++, ++k) {
		    colors->InsertComponent(k, 0, z(i, j));
		}
	    }
	    grid->GetPointData()->SetScalars(colors);

	    points->Delete();
	    colors->Delete();

	    return scale;

	}



	void renderer(vtkStructuredGrid *grid, double scale, vtkRenderer *rend,
		      bool draw_axes, bool draw_colorbar, bool draw_box,
		      bool do_warp, const double *range = NULL) {
	    //assert(has_data);

	    // filter to geometry primitive
	    vtkStructuredGridGeometryFilter *geometry =
		vtkStructuredGridGeometryFilter::New();
	    geometry->SetInput(grid);

	    // warp to fit in box
	    vtkWarpScalar *warp = vtkWarpScalar::New();
	    if (do_warp) {
		//double scale = Lxy/Lz;
		warp->SetInputConnection(geometry->GetOutputPort());
		warp->XYPlaneOn();
		warp->SetScaleFactor(scale);
	    }

	    // map gridfunction
	    vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
	    if (do_warp)
		mapper->SetInputConnection(warp->GetOutputPort());
	    else
		mapper->SetInputConnection(geometry->GetOutputPort());

	    double range_[2];
	    if (!range)  {
		grid->GetScalarRange(range_);
		range = range_;
	    }
	    mapper->SetScalarRange(range[0], range[1]);//tmp[0], tmp[1]);

	    // create plot surface actor
	    vtkActor *surfplot = vtkActor::New();
	    surfplot->SetMapper(mapper);

	    // create outline
	    vtkOutlineFilter *outlinefilter = vtkOutlineFilter::New();
	    if (do_warp)
		outlinefilter->SetInputConnection(warp->GetOutputPort());
	    else
		outlinefilter->SetInputConnection(geometry->GetOutputPort());
	    vtkPolyDataMapper *outlineMapper = vtkPolyDataMapper::New();
	    outlineMapper->SetInput(outlinefilter->GetOutput());
	    vtkActor *outline = vtkActor::New();
	    outline->SetMapper(outlineMapper);
	    outline->GetProperty()->SetColor(0, 0, 0);

	    // create axes
	    vtkAxesActor* axes = vtkAxesActor::New();
	    axes->SetShaftTypeToCylinder();
	    axes->SetNormalizedShaftLength( 0.85, 0.85, 0.85);
	    axes->SetNormalizedTipLength( 0.15, 0.15, 0.15);
	    axes->SetCylinderRadius( 0.500 * axes->GetCylinderRadius() );
	    axes->SetConeRadius( 1.025 * axes->GetConeRadius() );
	    axes->SetSphereRadius( 1.500 * axes->GetSphereRadius() );
	    vtkTextProperty* text_prop_ax = axes->GetXAxisCaptionActor2D()->
		GetCaptionTextProperty();
	    text_prop_ax->SetColor(0.0, 0.0, 0.0);
	    text_prop_ax->SetFontFamilyToArial();
	    text_prop_ax->SetFontSize(8);
	    axes->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->
		ShallowCopy(text_prop_ax);
	    axes->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->
		ShallowCopy(text_prop_ax);

	    // create colorbar
	    vtkScalarBarActor *colorbar = vtkScalarBarActor::New();
	    colorbar->SetLookupTable(mapper->GetLookupTable());
	    colorbar->SetWidth(0.085);
	    colorbar->SetHeight(0.9);
	    colorbar->SetPosition(0.9, 0.1);
	    vtkTextProperty* text_prop_cb = colorbar->GetLabelTextProperty();
	    text_prop_cb->SetColor(1.0, 1.0, 1.0);
	    colorbar->SetLabelTextProperty(text_prop_cb);

	    // renderer
	    //rend->Clear();
	    rend->RemoveAllViewProps();

	    rend->AddActor(surfplot);
	    if (draw_box)
		rend->AddActor(outline);
	    if (draw_axes)
		rend->AddActor(axes);
	    if (draw_colorbar)
		rend->AddActor(colorbar);

	    rend->SetBackground(0.25, 0.25, 0.25);

	    // renderer is now set up!

	    // clean up
	    colorbar->Delete();
	    warp->Delete();
	    axes->Delete();
	    outline->Delete();
	    outlinefilter->Delete();
	    outlineMapper->Delete();
	    surfplot->Delete();
	    mapper->Delete();
	    geometry->Delete();

	}


    }
}

using namespace matrix::plot;

void Surface::initialize() {
    //std::cout << "start\n";
    renderer_ = new Renderer(800, 600);
}

void Surface::render_(const expression &E) {
    index_vector x(E.size1()), y(E.size2());
    const double *range = (range_[0] - range_[1]) ? range_ : NULL;
    double scale = plot::geometry(x, y, E, *renderer_);
    plot::renderer(*renderer_, scale, *renderer_, false, true, true, false, range);
    renderer_->render();
    //if (!renderer_->started()) thread_ = boost::thread(&Renderer::start, renderer_);
}

void matrix::plot::surface_(const expression &E) {
    index_vector x(E.size1()), y(E.size2());
    matplot::Surf_VTK p;
    p.surf(x,y,E);
}
