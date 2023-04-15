#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <memory>
#include <utility>
#include <algorithm>
#include <map>
#include <string>

#define STB_IMAGE_IMPLEMENTATION
#include "include/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "include/stb_image_write.h"

//______________________________________________________IMAGE___________________________________________________

template <typename pixel_type>
class Image 
{
	public:
		const size_t height_, width_, channels_;
	
	private:
		size_t size_;
		std::vector<pixel_type> pixel_data_;
	
	public:
		Image(size_t height, size_t width, size_t channels):
				
				height_(height),
				width_(width),
				channels_(channels),
				size_(height_ * width_ * channels_),
				pixel_data_(size_, pixel_type()) {}

		pixel_type& operator()(size_t x, size_t y, size_t z) 
		{
			size_t idx = x * width_ * channels_ + y * channels_ + z;

			return pixel_data_[idx];
		}

    friend void Save8bitRgbImage(const std::string&, const Image<uint8_t>&);
};

void Save8bitRgbImage(const std::string& path, const Image<uint8_t>& image) 
{
    stbi_write_png(path.c_str(), image.width_, image.height_, 
				STBI_rgb, image.pixel_data_.data(), 
				image.width_ * image.channels_);
}
//_____________________________________________SUPPORT_FUNC_DISTANCE____________________________________________
struct RGBColor 
{
    uint8_t r, g, b;
};

RGBColor MixColors(RGBColor first, RGBColor second, double alpha) 
{
    return 
    {
		uint8_t(sqrt(double(first.r) * double(first.r) * alpha + double(second.r) * double(second.r) * (1 - alpha))),
		uint8_t(sqrt(double(first.g) * double(first.g) * alpha + double(second.g) * double(second.g) * (1 - alpha))),
		uint8_t(sqrt(double(first.b) * double(first.b) * alpha + double(second.b) * double(second.b) * (1 - alpha)))
    };
}

class Color 
{
    RGBColor base_;
    RGBColor gradient_to_;
    RGBColor border_;
    double frequency_;
    double thickness_;

    bool has_gradient_;
    bool has_border_;

public:
    Color(RGBColor base):
			
			base_(base),
			has_gradient_(false),
			has_border_(false){}

    Color(RGBColor base, RGBColor gradient_to, double frequency = 1.0):
			
			base_(base),
			gradient_to_(gradient_to),
			frequency_(frequency),
			has_gradient_(true),
			has_border_(false){}

    Color(RGBColor base, double thickness, RGBColor border):
			
				base_(base),
				border_(border),
				thickness_(thickness),
				has_gradient_(false),
				has_border_(true){}

    Color(RGBColor base, RGBColor gradient_to, RGBColor border, double frequency=1.0, double thickness=0.1):
				
				base_(base),
				gradient_to_(gradient_to),
				border_(border),
				frequency_(frequency),
				thickness_(thickness),
				has_gradient_(true),
				has_border_(true){}

    RGBColor getColor(double distance=0.0, double arg=0.0) 
    {
        if(has_border_ && std::abs(distance) < thickness_) 
        {
            return border_;
        } 
        
        else if(has_gradient_) 
        {
            return MixColors(base_, gradient_to_, std::pow(std::sin(arg * frequency_), 2.0));
        }
        
        else 
        {
            return base_;
        }
    }
};

class SDF 
{
	public:
		SDF() = default;
		virtual double distance(double x, double y) = 0;
		virtual ~SDF() = default;
		virtual RGBColor getColor(double x, double y) = 0;
};

class Circle: public SDF 
{
    double x_, y_;
    double radius_;
    Color color_;

	public:
		Circle(double x, double y, double radius, Color color): x_(x), y_(y), radius_(radius), color_(color) {}
		double distance(double x, double y) override 
		{
			return std::sqrt((x - x_) * (x - x_) + (y - y_) * (y - y_)) - radius_;
		}
    
		RGBColor getColor(double x, double y) override 
		{
			return color_.getColor(distance(x, y), x - x_ - radius_);
		}
};

class AxisAlignedRectangle: public SDF 
{
    double x_, y_;
    double width_, height_;
    Color color_;

	public:
		AxisAlignedRectangle(double x, double y, double width, double height, Color color):
			
				x_(x),
				y_(y),
				width_(width),
				height_(height),
				color_(color){}

			double distance(double x, double y) override 
			{
				double dx = std::abs(x - x_) - width_;
				double dy = std::abs(y - y_) - height_;

				return std::sqrt(std::max(dx, 0.0) * std::max(dx, 0.0) 
				+ std::max(dy, 0.0) * std::max(dy, 0.0)) 
				+ std::min(std::max(dx, dy), 0.0);
			}

    RGBColor getColor(double x, double y) override 
    {
        return color_.getColor(distance(x, y), y - y_ - height_);
    }
};

class Segment: public SDF 
{
    double a_x_, a_y_, b_x_, b_y_;
    RGBColor color_;
	
	public:
		Segment(double a_x, double a_y, double b_x, double b_y, RGBColor color):
			
			a_x_(a_x),
			a_y_(a_y),
			b_x_(b_x),
			b_y_(b_y),
			color_(color){}

		double distance(double x, double y) override 	
		{
			double dx = x - a_x_;
			double dy = y - a_y_;
			double bax = b_x_ - a_x_;
			double bay = b_y_ - a_y_;

			double h = std::clamp((dx * bax + dy * bay) / (bax * bax + bay * bay), 0.0, 1.0);
			return std::sqrt((dx - bax * h) * (dx - bax * h) + (dy - bay * h) * (dy - bay * h));
		}

		RGBColor getColor(double x, double y) override 
		{
			return color_;
		}
};

class AxisAlignedEquilateralTriangle: public SDF 
{
    double x_, y_;
    double radius_;
    Color color_;

	public:
		AxisAlignedEquilateralTriangle(double x, double y, double radius, Color color):
				
				x_(x),
				y_(y),
				radius_(radius * 2 / std::sqrt(3)),
				color_(color){}

		double distance(double x, double y) override 
		{	
			double k = std::sqrt(3.0);
			
			double dx = std::abs(x - x_) - radius_;
			double dy = y_ - y + radius_ / k + radius_ / std::sqrt(3) / 2;

			if(dx + k * dy > 0.0) 
			{
				double tmp = dx;
				dx = (dx - k * dy) / 2.0;
				dy = (-k * tmp - dy) / 2.0;
			}
			
			dx -= std::clamp(dx, -2.0 * radius_, 0.0);

			return (dy > 0.0 ? -1.0 : 1.0) * std::sqrt(dx * dx + dy * dy);
		}

		RGBColor getColor(double x, double y) override 
		{
			return color_.getColor(distance(x, y), y - y_ - radius_);
		}
};

class SDFImage: public SDF 
{
    double x_, y_;
    double scale_;
    Color color_;

    int width_, height_;
    int max_side_;
    uint8_t* data_;

    uint8_t getPixel(size_t i, size_t j) 
    {
        if(i < height_ && j < width_) 
        {
            return data_[i * width_ + j];
        } 
        
        else 
        {
            return 0;
        }
    }

	public:
		
		SDFImage(const std::string& filepath, 
		double x, double y, double scale, 
		Color color): 
			x_(x), y_(y), scale_(scale),
			color_(color) 
			{
				int channels;
				data_ = stbi_load(filepath.c_str(), &width_, &height_, &channels, 1);
				max_side_ = std::max(width_, height_);
				
				if(data_ == nullptr) 
				{
					std::cerr << "Изображение не удалось загрузить" << std::endl;
					exit(1);
				}
				
				x_ = x_ - scale_ * width_ / max_side_ / 2;
				y_ = y_ - scale_ * height_ / max_side_ / 2;
			}

		double distance(double x, double y) override 
		{
			double ix = x - x_;
			double iy = y - y_;
			if(ix < 0 || iy < 0 || ix >= scale_ * width_ / max_side_ 
			|| iy >= scale_ * height_ / max_side_) 
			{
				return 1.0;
			} 
			
			else 
			{
            
				ix *= max_side_ / scale_;
				iy *= max_side_ / scale_;
				double interp_x = std::modf(ix, &ix);
				double interp_y = std::modf(iy, &iy);

				double pixels[4] = 
				{
                    static_cast<double>(getPixel(static_cast<size_t>(iy    ), static_cast<size_t>(ix    ))),
                    static_cast<double>(getPixel(static_cast<size_t>(iy    ), static_cast<size_t>(ix + 1))),
                    static_cast<double>(getPixel(static_cast<size_t>(iy + 1), static_cast<size_t>(ix    ))),
                    static_cast<double>(getPixel(static_cast<size_t>(iy + 1), static_cast<size_t>(ix + 1)))
				};
				double value = pixels[0] * (1 - interp_y) * (1 - interp_x) +
							pixels[1] * (1 - interp_y) * (    interp_x) +
							pixels[2] * (    interp_y) * (1 - interp_x) +
							pixels[3] * (    interp_y) * (    interp_x);
				
				return (128. - value) / 255.;
			}
		}

		RGBColor getColor(double x, double y) override 
		{
			return color_.getColor(distance(x, y), y - y_ - scale_);
		}
};

double sminCubic(double a, double b, double k)
{
    double h = std::max( k-abs(a-b), 0.0 ) / k;
    double m = h * h * h * 0.5;
    double s = m * k * (1.0/3.0);
    return (a<b) ? a-s : b-s;
}

double sminCubicCol(double a, double b, double k)
{
    double h = std::max( k - abs(a-b), 0.0 )/k;
    double m = h * h * h * 0.5;
    double s = m * k * (1.0/3.0);
    return (a<b) ? m : 1-m;
}

class Intersection: public SDF 
{
    std::shared_ptr<SDF> first_;
    std::shared_ptr<SDF> second_;
    bool smooth_;
    double smoothness_;

	public:
		Intersection(std::shared_ptr<SDF>&& first, std::shared_ptr<SDF>&& second,
					bool smooth=false, double smoothness=0.125):
         
			first_(std::move(first)),
			second_(std::move(second)),
			smooth_(smooth),
			smoothness_(smoothness){}

		double distance(double x, double y) override 
		{
			if(smooth_) 
			{
				return sminCubic(first_->distance(x, y), second_->distance(x, y), smoothness_);
			} 
			
			else 
			{
				return std::min(first_->distance(x, y), second_->distance(x, y));
			}
		}

		RGBColor getColor(double x, double y) override 
		{
			double first_dist = first_->distance(x, y);
			double second_dist = second_->distance(x, y);
        
			if(smooth_) 
			{
				double blend = sminCubicCol(second_dist, first_dist, smoothness_);
				auto fcol = first_->getColor(x, y);
				auto scol = second_->getColor(x, y);
				return MixColors(fcol, scol, blend);
			} 
        
			else 
			{
				if(first_dist < second_dist) 
				{
					return first_->getColor(x, y);
				} 
				
				else 
				{
					return second_->getColor(x, y);
				}
			}
		}
};

class Overlay: public SDF 
{
    std::shared_ptr<SDF> top_;
    std::shared_ptr<SDF> bottom_;
    double alpha_;

	public:
		Overlay(std::shared_ptr<SDF>&& top, std::shared_ptr<SDF>&& bottom, double alpha=0.5):
			
			top_(std::move(top)),
			bottom_(std::move(bottom)),
			alpha_(alpha){}

		double distance(double x, double y) override 
		{
			return std::min(top_->distance(x, y), bottom_->distance(x, y));
		}

		RGBColor getColor(double x, double y) override 
		{
			double top_dist = top_->distance(x, y);
			double bottom_dist = bottom_->distance(x, y);

			RGBColor top_color = top_->getColor(x, y);
			RGBColor bottom_color = bottom_->getColor(x, y);

			double eps = 2e-3;
			if(top_dist < eps && bottom_dist < eps) 
			{
				return MixColors(top_color, bottom_color, alpha_);
			} 
			
			else if(bottom_dist < eps) 
			{
				return bottom_color;
			} 
			
			else 
			{
				return top_color;
			}
		}
};
//_______________________________________________SCENE_________________________________________________

class Scene 
{
    std::vector<std::shared_ptr<SDF>> objects_;             
    double x_min_, x_max_, y_min_, y_max_; 
    RGBColor background_;
	
	public:
		Scene(const std::vector<std::shared_ptr<SDF>>& objects, 
		double x_min, double x_max, double y_min, double y_max, RGBColor background):
					
					objects_(objects),
					x_min_(x_min),
					x_max_(x_max),
					y_min_(y_min),
					y_max_(y_max),
					background_(background){}

    template<typename pixel_type>
    void RenderToImage(Image<pixel_type>& image, double eps=1e-3) 
    {
        for(int i = 0; i < image.height_; ++i) 
        {
            for(int j = 0; j < image.width_; ++j) 
            {
                double y = y_min_ + double(i) / image.width_ * (y_max_ - y_min_);
                double x = x_min_ + double(j) / image.height_ * (x_max_ - x_min_);

                std::shared_ptr<SDF> hit;
                for(const auto& object : objects_) 
                {
                    if(object->distance(x, y) < eps) 
                    {
                        hit = object;
                        break;
                    }
                }
                
                RGBColor color = background_;
                
                if(hit) 
                {
                    color = hit->getColor(x, y);
                }
                
                image(i, j, 0) = color.r;
                image(i, j, 1) = color.g;
                image(i, j, 2) = color.b;
            }
        }
    }
};

//________________________________________________HEAD__________________________________________________
//-----------------------1------------------------
Scene Scene1() 
{
	//x y sizer
	//red green blue
    return Scene(
    {
		
        std::make_shared<Circle>(0.66, 0.0, 0.25, Color({225, 11, 16})),
        std::make_shared<Circle>(-0.75, -0.75, 0.15, Color({0, 255, 16})),
        
        std::make_shared<AxisAlignedRectangle>(-0.66, 0.0, 0.25, 0.25, Color({225, 11, 16})),
        
        std::make_shared<AxisAlignedRectangle>(0.75, 0.75, 0.15, 0.15, Color({90, 180, 200})),
        std::make_shared<AxisAlignedRectangle>(0.80, 0.60, 0.40, 0.15, Color({255, 20, 180})),
        std::make_shared<AxisAlignedEquilateralTriangle>(0.0, -0.30, 0.25, Color({255, 182, 70})),
        
        std::make_shared<AxisAlignedRectangle>(0.66, 0.0, 0.25, 0.25, Color({0, 0, 0})),}, 
        -1, 1, -1, 1, {192, 192, 192});
}

//-----------------------2------------------------
Scene Scene2() 	
{
    return Scene(
    {
		std::make_shared<Overlay>(
		std::make_shared<SDFImage>("../aaa.png", -0.7, -0.7, 1., Color({255, 190, 180})),
		std::make_shared<SDFImage>("../ppp.png", -0.625, -0.6, 1., Color({255, 255, 255})),0.7),
		
		std::make_shared<SDFImage>("../amogus.png", 0.7, 0.7, 0.3, Color({50, 20, 90})),
        
        std::make_shared<SDFImage>("../sdf.png",-0.25, 0.1, 1., Color({22, 190, 255}, {0, 0, 255}, {255, 182, 70}, 3, 0.4)),
        std::make_shared<SDFImage>("../sdf.png",0.25, -0.25, 0.30, Color({120, 180, 90}, {0, 255, 190}, {255, 182, 70}, 3, 0.4)),
        std::make_shared<Intersection>(std::make_shared<Circle>(0.725, 0.7, 0.125, Color({0, 255, 0})),
		std::make_shared<Intersection>(std::make_shared<Circle>(0.475, 0.7, 0.125, Color({255, 0, 0})),
		
		std::make_shared<Circle>(0.6, 0.5, 0.125, Color({0, 0, 255})),true),true)}, 
		-1, 1, -1, 1, {192, 192, 192});
}
//-----------------------------3-----------------------
Scene Scene3()
{
	return Scene(
	{
			std::make_shared<Overlay>(
            std::make_shared<AxisAlignedRectangle>(0., 0., 0.40, 0.40, Color({255, 80, 100}, {100, 255, 255}, {100, 30, 200}, 4, 0.03)),
            std::make_shared<Intersection>(std::make_shared<Intersection>(std::make_shared<Intersection>(
			std::make_shared<Intersection>(std::make_shared<Intersection>(std::make_shared<Intersection>(
			std::make_shared<Intersection>(std::make_shared<Intersection>(std::make_shared<Intersection>(std::make_shared<Intersection>(
			std::make_shared<AxisAlignedRectangle>(0., 0., 0.55, 0.75, Color({120, 200, 200}, {90, 120, 255}, {255,255,255}, 2, 0.1)),
			
			std::make_shared<Segment>(-0.75, -0.7, 0.75, -0.45, RGBColor{120, 0, 200}), true),
			std::make_shared<Segment>(0.75, 0.7, 0.75, -0.5, RGBColor{255, 120, 120}), true),
			std::make_shared<Segment>(0.50, 0.9, 0.75, -0.5, RGBColor{150, 220, 200}), true),
			std::make_shared<Segment>(-0.75, -0.7, -0.75, -0.5, RGBColor{100, 100, 100}), true),
			std::make_shared<Segment>(-0.50, -0.9, -0.75, -0.5, RGBColor{200, 90, 40}), true),
			std::make_shared<Segment>(-0.75, -0.1, 0.75, -0.25, RGBColor{42, 42, 42}), true),
			std::make_shared<Segment>(-0.75, 0.1, 0.75, 0.3, RGBColor{255, 10, 120}), true),
			std::make_shared<Segment>(0.75, -0.7, -0.75, -0.5, RGBColor{80, 80, 80}), true),
			std::make_shared<Segment>(0.50, -0.9, -0.75, -0.5, RGBColor{110, 88, 130}), true),
			std::make_shared<Segment>(-0.75, 0.7, 0.75, 0.45, RGBColor{228, 88, 14}),true))
    }
    , -1, 1, -1, 1, {192, 192, 192});
}

int main() 
{
    Image<uint8_t> rgb_image(1024, 1024, 3);
	//---------------------------------------------1-----------------------------------------------------------------------
    auto scene = Scene1();
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    scene.RenderToImage(rgb_image, 2e-3);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    Save8bitRgbImage("../scene1.png", rgb_image);
    std::cout << "Scene 1 Times: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() 
    << "ms" << std::endl;

//---------------------------------------------2-----------------------------------------------------------------------
    scene = Scene2();
    
    begin = std::chrono::steady_clock::now();
    scene.RenderToImage(rgb_image, 2e-3);
    end = std::chrono::steady_clock::now();
    
    Save8bitRgbImage("../scene2.png", rgb_image);
    std::cout << "Scene 2 Times: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() 
    << "ms" << std::endl;
//---------------------------------------------2-----------------------------------------------------------------------

    scene = Scene3();
    
    begin = std::chrono::steady_clock::now();
    scene.RenderToImage(rgb_image, 2e-3);
    end = std::chrono::steady_clock::now();
    
    Save8bitRgbImage("../scene3.png", rgb_image);
   
    std::cout << "Scene 3 Times: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() 
    << "ms" << std::endl;
}
