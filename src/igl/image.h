#ifndef _IMAGE_H_
#define _IMAGE_H_

#include "common/common.h"
#include "vmath/vmath.h"

///@file igl/image.h Images. @ingroup igl
///@defgroup image Images
///@ingroup igl
///@{

/// generic image
template<typename T>
struct image {
    /// Default Constructor (empty image)
	image() : _w(0), _h(0) { }
    /// Size Constructor (sets width and height)
	image(int w, int h) : _w(w), _h(h), _d(w*h,T()) { }
    
    /// image width
    int width() const { return _w; }
    /// image height
    int height() const { return _h; }
	
    /// element access
	T& at(int i, int j) { return _d[j*_w+i]; }
    /// element access
	const T& at(int i, int j) const { return _d[j*_w+i]; }
    
    /// data access
    T* data() { return _d.data(); }
    /// data access
    const T* data() const { return _d.data(); }

    /// set image elements
    void set(const T& v) { for(auto& d : _d) d = v; }

private:
	int _w, _h;
	vector<T> _d;
};

///@name image typedefs
///@{
typedef image<vec3f> image3f;
typedef image<vec3ub> image3ub;
///@}

///@name image operations
///@{
template<typename T>
inline void copy_image(image<T>& dest, image<T>& source) {
    error_if_not(dest.with() == source.width() and dest.height() == source.height(), "images should be same size");
    for(int y = 0; y < dest.height(); y ++) {
        for(int x = 0; x < dest.width(); x ++) {
            dest.at(x,y) = source.at(x,y);
        }
    }
}

template<typename T>
inline void copy_subimage(image<T>& dest, image<T>& source, int ox, int oy) {
    for(int y = 0; y < dest.height(); y ++) {
        for(int x = 0; x < dest.width(); x ++) {
            dest.at(x+ox,y+oy) = source.at(x,y);
        }
    }
}

template <typename T1, typename T2>
inline image<T1> image_cast(const image<T2>& img) {
    image<T1> ret(img.width(), img.height());
    for(int j = 0; j < img.height(); j ++) {
        for(int i = 0; i < img.width(); i ++) {
            
        }
    }
    return ret;
}

template <typename T>
inline image<T> flipy(const image<T>& img) {
    image<T> ret(img.width(),img.height());
    for(int j = 0; j < img.height(); j ++) {
        for(int i = 0; i < img.width(); i ++) {
            ret.at(i,img.height()-1-j) = img.at(i,j);
        }
    }
    return ret;
}
///@}

///@name image io
///@{
void imageio_write_pfm(const string& filename, const image<vec3f>& img, bool flipY = false);
void imageio_write_ppm(const string& filename, const image<vec3f>& img, bool flipY = false);
void imageio_write_png(const string& filename, const image<vec3f>& img, bool flipY = false);

image<vec3f> imageio_read_pnm3f(const string& filename);
image<vec3f> imageio_read_png3f(const string& filename);
///@}

///@}

#endif
