#include "tesselate.h"

///@file igl/tesselate.cpp Tesselation. @ingroup igl

/// Tesselate parametric Shape (two parameters)
/// @param shape_frame Function that returns the frame of the Shape at the given vec2f coordinate
/// @param ur Number of divisions along first parameter
/// @param vr Number of divisions along second parameter
/// @param ccw Specifies whether faces are to be generated counter-clockwise or clockwise
/// @param smooth Specifies whether the tesselated surface is to have smooth normals
/// @return The tesselated shape
Shape* _tesselate_shape_uniform(const function<frame3f (const vec2f&)> shape_frame,
                                int ur, int vr, bool ccw, bool smooth) {
    auto tesselation = new Mesh();
    auto line = vector<vec2i>();
    
    for(int i = 0; i <= ur; i ++) {
        for(int j = 0; j <= vr; j ++) {
            float u = i / float(ur);
            float v = j / float(vr);
            vec2f uv = vec2f(u,v);
            auto f = shape_frame(uv);
            tesselation->pos.push_back(f.o);
            if(smooth) tesselation->norm.push_back(f.z);
        }
    }
    
    for(int i = 0; i < ur; i ++) {
        for(int j = 0; j < vr; j ++) {
            vec4i f = vec4i((i+0)*(vr+1)+(j+0),(i+0)*(vr+1)+(j+1),(i+1)*(vr+1)+(j+1),(i+1)*(vr+1)+(j+0));
            if(not ccw) { swap(f.y,f.w); }
            tesselation->quad.push_back(f);
        }
    }
    
    return tesselation;
}

/// Tesselate Shape (uniformly) with lines
/// @param shape_frame Function that returns the frame of the Shape at the given vec2f coordinate
/// @param ur Number of divisions along first parameter for Shape
/// @param vr Number of divisions along second parameter for shape
/// @param ul Number of divisions along first parameter for lines
/// @param vl Number of divisions along second parameter for lines
/// @param ccw Specifies whether faces are to be generated counter-clockwise or clockwise
/// @param smooth Specifies whether the tesselated surface is to have smooth normals
/// @return The tesselated shape with lines
Shape* _tesselate_shape_uniform(const function<frame3f (const vec2f&)> shape_frame,
                                int ur, int vr, int ul, int vl, bool ccw, bool smooth) {
    auto tesselation = cast<Mesh>(_tesselate_shape_uniform(shape_frame, ur, vr, ccw, smooth));
    auto line = vector<vec2i>();
    
    for(int li = 0; li <= ul; li ++) {
        int i = li * ur / ul;
        for(int j = 0; j < vr; j ++) {
            vec2i l = vec2i((i+0)*(vr+1)+(j+0),(i+0)*(vr+1)+(j+1));
            tesselation->_tesselation_lines.push_back(l);
        }
    }
    for(int lj = 0; lj <= vl; lj ++) {
        int j = lj * vr / vl;
        for(int i = 0; i < ur; i ++) {
            vec2i l = vec2i((i+0)*(vr+1)+(j+0),(i+1)*(vr+1)+(j+0));
            tesselation->_tesselation_lines.push_back(l);
        }
    }
    
    return tesselation;
}

/// Tesselate parametric Shape (one parameter) that has a radius.  ex: bezier spline
/// @param shape_frame Function that returns the frame of Shape at the given single float coordinate
/// @param shape_radius Function that returns the radius of Shape at the given coordinate
/// @param ur Number of divisions along parameter
/// @param smooth Specifies whether the tesselated surface is to have smooth normals
/// @return The tesselated Shape
Shape* _tesselate_shape_uniform(const function<frame3f (float)> shape_frame,
                                const function<float (float)> shape_radius,
                                int ur, bool smooth) {
    auto tesselation = new LineSet();

    put_your_code_here("Bezier Spline Tesselation");
    
    return tesselation;
}

LineSet* _tesselate_lineset_once(LineSet* lines) {
    auto tesselation = new LineSet();
    
    // add vertices
    tesselation->pos = lines->pos;
    tesselation->radius = lines->radius;
    
    // add edge vertices
    int lvo = tesselation->pos.size();
    for(auto l : lines->line) {
        tesselation->pos.push_back((lines->pos[l.x]+lines->pos[l.y])/2);
        tesselation->radius.push_back((lines->radius[l.x]+lines->radius[l.y])/2);
    }
    
    // add lines
    for(int lid = 0; lid < lines->line.size(); lid ++) {
        auto l = lines->line[lid];
        tesselation->line.push_back(vec2i(l.x,lid+lvo));
        tesselation->line.push_back(vec2i(lid+lvo,l.y));
    }
    
    return tesselation;
}

TriangleMesh* _tesselate_trianglemesh_once(TriangleMesh* mesh) {
    auto tesselation = new TriangleMesh();
    
    // adjacency
    auto adj = EdgeHashTable(mesh->triangle,vector<vec4i>());
    
    // add vertices
    tesselation->pos = mesh->pos;
    tesselation->norm = mesh->norm;
    tesselation->texcoord = mesh->texcoord;
    
    // add edge vertices
    int evo = tesselation->pos.size();
    for(auto e : adj.edges) {
        tesselation->pos.push_back(mesh->pos[e.x]*0.5+mesh->pos[e.y]*0.5);
        if(not tesselation->norm.empty()) tesselation->norm.push_back(normalize(mesh->norm[e.x]*0.5+mesh->norm[e.y]*0.5));
        if(not tesselation->texcoord.empty()) tesselation->texcoord.push_back(mesh->texcoord[e.x]*0.5+mesh->texcoord[e.y]*0.5);
    }
    
    // add triangles
    for(auto f : mesh->triangle) {
        auto ve = vec3i(adj.edge(f.x, f.y),adj.edge(f.y, f.z),adj.edge(f.z, f.x))+vec3i(evo,evo,evo);
        tesselation->triangle.push_back(vec3i(f.x,ve.x,ve.z));
        tesselation->triangle.push_back(vec3i(f.y,ve.y,ve.x));
        tesselation->triangle.push_back(vec3i(f.z,ve.z,ve.y));
        tesselation->triangle.push_back(ve);
    }
    
    // add lines
    if(mesh->_tesselation_lines.size() > 0) {
        for(auto l : mesh->_tesselation_lines) {
            int ve = adj.edge(l.x, l.y)+evo;
            tesselation->_tesselation_lines.push_back(vec2i(l.x,ve));
            tesselation->_tesselation_lines.push_back(vec2i(ve,l.y));
        }
    }
    
    return tesselation;
}

/// Tesselate Mesh one time
/// @param mesh The Mesh to tesselate
/// @return The tesselated Mesh
Mesh* _tesselate_mesh_once(Mesh* mesh) {
    auto tesselation = new Mesh();
    
    // adjacency
    auto adj = EdgeHashTable(mesh->triangle, mesh->quad);
    
    // add vertices
    tesselation->pos = mesh->pos;

    // Apparently, the below aren't needed
    //tesselation->norm = mesh->norm;
    //tesselation->texcoord = mesh->texcoord;
    
    // add edge vertices
    int evo = tesselation->pos.size();
    for(auto e : adj.edges) {
        tesselation->pos.push_back(mesh->pos[e.x]*0.5+mesh->pos[e.y]*0.5);
        // 80char limit yall
        if (not tesselation->norm.empty()) {
            tesselation->norm.push_back(
                normalize(mesh->norm[e.x]*0.5+mesh->norm[e.y]*0.5)

            );
        }
        if (not tesselation->texcoord.empty()) {
            tesselation->texcoord.push_back(
                mesh->texcoord[e.x]*0.5+mesh->texcoord[e.y]*0.5

            );
        }
    }
    
    // add triangles
    for(auto f : mesh->triangle) {
        auto ve = vec3i(
                    adj.edge(f.x, f.y),
                    adj.edge(f.y, f.z),
                    adj.edge(f.z, f.x)) + vec3i(evo, evo, evo);
        tesselation->triangle.push_back(vec3i(f.x, ve.x, ve.z));
        tesselation->triangle.push_back(vec3i(f.y, ve.y, ve.x));
        tesselation->triangle.push_back(vec3i(f.z, ve.z, ve.y));
        tesselation->triangle.push_back(ve);
    }
    
    // add quads
    for(auto f : mesh->quad) {
        auto ve = vec4i(
                    adj.edge(f.x, f.y),
                    adj.edge(f.y, f.z),
                    adj.edge(f.z, f.w),
                    adj.edge(f.w, f.x)) + vec4i(evo, evo, evo, evo);

        // Create midpoint geometry as centroid of vertices
        tesselation->pos.push_back( 
            (mesh->pos[f.x]+mesh->pos[f.y]+mesh->pos[f.z]+mesh->pos[f.w]) / 4

        );

        // Index of newly created vertex
        auto cent_ind = tesselation->pos.size() - 1;

        tesselation->quad.push_back(vec4i(f.x, ve.x, cent_ind, ve.w));
        tesselation->quad.push_back(vec4i(f.y, ve.y, cent_ind, ve.x));
        tesselation->quad.push_back(vec4i(f.z, ve.z, cent_ind, ve.y));
        tesselation->quad.push_back(vec4i(f.w, ve.w, cent_ind, ve.z));

    }

    // add lines
    if(mesh->_tesselation_lines.size() > 0) {
        for(auto l : mesh->_tesselation_lines) {
            int ve = adj.edge(l.x, l.y)+evo;
            tesselation->_tesselation_lines.push_back(vec2i(l.x,ve));
            tesselation->_tesselation_lines.push_back(vec2i(ve,l.y));
        }
    }    
    return tesselation;

}

/// Apply Catmull-Clark subdivision surface rules on subdiv
/// @param subdiv The mesh to subdivide
/// @return The subdivided mesh
CatmullClarkSubdiv* _tesselate_catmullclark_once(CatmullClarkSubdiv* subdiv) {
    int i;

    message("Running _tesselate_catmullclarck_once(), which is currently incomplete");

    // Part 1: Basically a direct copy of _tesselate_mesh_once(), constructs
    //     tesselations for entire subdiv mesh
    auto tesselation = new CatmullClarkSubdiv();    

    // adjacency
    // No triangles in a subdivm, only quadrilaterals
    auto adj = EdgeHashTable(vector<vec3i>(), subdiv->quad);
    
    // add vertices
    tesselation->pos = subdiv->pos;

    // Apparently, the below aren't needed
    //tesselation->norm = subdiv->norm;
    //tesselation->texcoord = subdiv->texcoord;
    
    // add edge vertices
    int evo = tesselation->pos.size();
    for(auto e : adj.edges) {
        tesselation->pos.push_back(subdiv->pos[e.x]*0.5+subdiv->pos[e.y]*0.5);
        // 80char limit yall
        if (not tesselation->norm.empty()) {
            tesselation->norm.push_back(
                normalize(subdiv->norm[e.x]*0.5+subdiv->norm[e.y]*0.5)

            );
        }
        if (not tesselation->texcoord.empty()) {
            tesselation->texcoord.push_back(
                subdiv->texcoord[e.x]*0.5+subdiv->texcoord[e.y]*0.5

            );
        }
    }
    
    // add quads
    for(auto f : subdiv->quad) {
        auto ve = vec4i(
                    adj.edge(f.x, f.y),
                    adj.edge(f.y, f.z),
                    adj.edge(f.z, f.w),
                    adj.edge(f.w, f.x)) + vec4i(evo, evo, evo, evo);

        // Create midpoint geometry as centroid of vertices
        tesselation->pos.push_back( 
            (subdiv->pos[f.x]+subdiv->pos[f.y]+subdiv->pos[f.z]+subdiv->pos[f.w]) / 4

        );

        // Index of newly created vertex
        auto cent_ind = tesselation->pos.size() - 1;

        tesselation->quad.push_back(vec4i(f.x, ve.x, cent_ind, ve.w));
        tesselation->quad.push_back(vec4i(f.y, ve.y, cent_ind, ve.x));
        tesselation->quad.push_back(vec4i(f.z, ve.z, cent_ind, ve.y));
        tesselation->quad.push_back(vec4i(f.w, ve.w, cent_ind, ve.z));

    }

    // add lines
    if(subdiv->_tesselation_lines.size() > 0) {
        for(auto l : subdiv->_tesselation_lines) {
            int ve = adj.edge(l.x, l.y)+evo;
            tesselation->_tesselation_lines.push_back(vec2i(l.x,ve));
            tesselation->_tesselation_lines.push_back(vec2i(ve,l.y));
        }
    }    
    
    // Part 2: Averaging step
    //     Adjusts existing vertices by the average of adjacent face centroids

    // Averaging holders, used in part 2;
    //     avg_v holds sum of vertices and avg_n holds number in sum
    vec3f avg_v[tesselation->pos.size()];
    int avg_n[tesselation->pos.size()];
    //printf("%d %d\n", tesselation->pos.size(), sizeof avg_v / sizeof avg_v[0]);
    //fflush(stdout);
    memset(avg_v, 0, sizeof avg_v / sizeof avg_v[0]);
    memset(avg_n, 0, sizeof avg_n / sizeof avg_n[0]);

    for (auto f : tesselation->quad) {
        // Find centroid of this quad face
        auto cent = (tesselation->pos[f.x] + tesselation->pos[f.y]
                    + tesselation->pos[f.z] + tesselation->pos[f.w])
                    / 4;

        // Go through all of the vertices in this face and average position
        avg_v[f.x] += cent;
        avg_n[f.x] += 1;
        avg_v[f.y] += cent;
        avg_n[f.y] += 1;
        avg_v[f.z] += cent;
        avg_n[f.z] += 1;
        avg_v[f.w] += cent;
        avg_n[f.w] += 1;

    }

    // Run the average calculation
    for (i = 0; i < tesselation->pos.size(); i++) {
        avg_v[i] /= avg_n[i];

    }

    // Part 3: Adjust by average
    //     
    for (i = 0; i < tesselation->pos.size(); i++) {
        avg_v[i] = tesselation->pos[i] + (avg_v[i] - tesselation->pos[i])
                   * (4/avg_n[i]);

    }
    tesselation->pos = vector<vec3f>(avg_v,
                                     avg_v + (sizeof avg_v / sizeof avg_v[0]));
    return tesselation;

}

/// Apply subdivision surface rules on subdiv
/// @param subdiv The mesh to subdivide
/// @return The subdivided mesh
Subdiv* _tesselate_subdiv_once(Subdiv* subdiv) {
    auto tesselation = new Subdiv();
    
    // linear subdivision like triangle mesh
    // adjacency
    auto adj = EdgeHashTable(subdiv->triangle,subdiv->quad);
    
    // add vertices
    tesselation->pos = subdiv->pos;
    tesselation->texcoord = subdiv->texcoord;
    
    // add edge vertices
    int evo = tesselation->pos.size();
    for(auto e : adj.edges) {
        tesselation->pos.push_back(subdiv->pos[e.x]*0.5+subdiv->pos[e.y]*0.5);
        if(not tesselation->texcoord.empty()) tesselation->texcoord.push_back(subdiv->texcoord[e.x]*0.5+subdiv->texcoord[e.y]*0.5);
    }
    
    // add face vertices
    int fvo = tesselation->pos.size();
    for(auto f : subdiv->quad) {
        tesselation->pos.push_back(subdiv->pos[f.x]*0.25+subdiv->pos[f.y]*0.25+subdiv->pos[f.z]*0.25+subdiv->pos[f.w]*0.25);
        if(not tesselation->texcoord.empty()) tesselation->texcoord.push_back(subdiv->texcoord[f.x]*0.25+subdiv->texcoord[f.y]*0.25+subdiv->texcoord[f.z]*0.25+subdiv->texcoord[f.w]*0.25);
    }
    
    // add triangles
    for(auto f : subdiv->triangle) {
        auto ve = vec3i(adj.edge(f.x, f.y),adj.edge(f.y, f.z),adj.edge(f.z, f.x))+vec3i(evo,evo,evo);
        tesselation->triangle.push_back(vec3i(f.x,ve.x,ve.z));
        tesselation->triangle.push_back(vec3i(f.y,ve.y,ve.x));
        tesselation->triangle.push_back(vec3i(f.z,ve.z,ve.y));
        tesselation->triangle.push_back(ve);
    }
    
    // add quads
    for(int fid = 0; fid < subdiv->quad.size(); fid ++) {
        auto f = subdiv->quad[fid];
        auto ve = vec4i(adj.edge(f.x, f.y),adj.edge(f.y, f.z),adj.edge(f.z, f.w),adj.edge(f.w, f.x))+vec4i(evo,evo,evo,evo);
        auto vf = fid+fvo;
        tesselation->quad.push_back(vec4i(f.x,ve.x,vf,ve.w));
        tesselation->quad.push_back(vec4i(f.y,ve.y,vf,ve.x));
        tesselation->quad.push_back(vec4i(f.z,ve.z,vf,ve.y));
        tesselation->quad.push_back(vec4i(f.w,ve.w,vf,ve.z));
    }
    
    // creases
    tesselation->crease_vertex = subdiv->crease_vertex;
    for(auto e : subdiv->crease_edge) {
        tesselation->crease_edge.push_back({e.x,evo+adj.edge(e.x, e.y)});
        tesselation->crease_edge.push_back({evo+adj.edge(e.x, e.y),e.y});
    }
    
    // add lines
    if(subdiv->_tesselation_lines.size() > 0) {
        for(auto l : subdiv->_tesselation_lines) {
            int ve = adj.edge(l.x, l.y)+evo;
            tesselation->_tesselation_lines.push_back(vec2i(l.x,ve));
            tesselation->_tesselation_lines.push_back(vec2i(ve,l.y));
        }
    }
    
    // mark creases
    auto cvertex = vector<bool>(tesselation->pos.size(),false);
    auto cedge = vector<bool>(tesselation->pos.size(),false);
    for(auto vid : tesselation->crease_vertex) cvertex[vid] = true;
    for(auto e : tesselation->crease_edge) for(auto vid : e) cedge[vid] = true;
    
    // averaging
    auto npos = vector<vec3f>(tesselation->pos.size(),zero3f);
    auto ntexcoord = vector<vec2f>(tesselation->texcoord.size(),zero2f);
    auto weight = vector<float>(tesselation->pos.size(),0);
    auto nquad = vector<int>(tesselation->pos.size(),0);
    auto ntriangle = vector<int>(tesselation->pos.size(),0);
    for(auto f : tesselation->quad) {
        for(auto vid : f) {
            if(cedge[vid] or cvertex[vid]) continue;
            auto w = pi/2;
            npos[vid] += (tesselation->pos[f.x]+tesselation->pos[f.y]+tesselation->pos[f.z]+tesselation->pos[f.w])*w/4;
            if(not ntexcoord.empty()) ntexcoord[vid] += (tesselation->texcoord[f.x]+tesselation->texcoord[f.y]+tesselation->texcoord[f.z]+tesselation->texcoord[f.w])*w/4;
            weight[vid] += w;
            nquad[vid] ++;
        }
    }
    for(auto f : tesselation->triangle) {
        for(int i = 0; i < 3; i ++) {
            auto vid = f[i]; auto vid1 = f[(i+1)%3]; auto vid2 = f[(i+2)%3];
            if(cedge[vid] or cvertex[vid]) continue;
            auto w = pi/3;
            npos[vid] += (tesselation->pos[vid]/4+tesselation->pos[vid1]*(3/8.0)+tesselation->pos[vid2]*(3/8.0))*w;
            if(not ntexcoord.empty()) ntexcoord[vid] += (tesselation->texcoord[vid]/4+tesselation->texcoord[vid1]*(3/8.0)+tesselation->texcoord[vid2]*(3/8.0))*w;
            weight[vid] += w;
            ntriangle[vid] ++;
        }
    }
    
    // handle creases
    for(auto e : tesselation->crease_edge) {
        for(auto vid : e) {
            if(cvertex[vid]) continue;
            npos[vid] += (tesselation->pos[e.x]+tesselation->pos[e.y])/2;
            if(not ntexcoord.empty()) ntexcoord[vid] += (tesselation->texcoord[vid]+tesselation->texcoord[vid])/2;
            weight[vid] += 1;
        }
    }
    for(auto vid : tesselation->crease_vertex) {
        npos[vid] += tesselation->pos[vid];
        if(not ntexcoord.empty()) ntexcoord[vid] += tesselation->texcoord[vid];
        weight[vid] += 1;
    }
    
    // normalization
    for(auto i : range(tesselation->pos.size())) {
        npos[i] /= weight[i];
        if(not ntexcoord.empty()) ntexcoord[i] /= weight[i];
    }
    
    // correction
    for(auto i : range(tesselation->pos.size())) {
        if(cedge[i] or cvertex[i]) continue;
        float w = 0;
        if(tesselation->quad.empty()) w = 5/3.0 - (8/3.0)*pow(3/8.0+1/4.0*cos(2*pi/ntriangle[i]),2);
        else if(tesselation->triangle.empty()) w = 4.0f / nquad[i];
        else w = (nquad[i] == 0 and ntriangle[i] == 3) ? 1.5f : 12.0f / (3 * nquad[i] + 2 * ntriangle[i]);
        npos[i] = tesselation->pos[i] + (npos[i] - tesselation->pos[i])*w;
        if(not ntexcoord.empty()) ntexcoord[i] = tesselation->texcoord[i] + (ntexcoord[i] - tesselation->texcoord[i])*w;
    }
    
    // set tesselation back
    tesselation->pos = npos;
    tesselation->texcoord = ntexcoord;
    
    return tesselation;
}

Shape* _tesselate_recursive(const function<Shape*(Shape*)>& tesselate_once, Shape* tesselation, int level, bool smooth,
                            const function<Shape*(Shape*)>& to_mesh = function<Shape*(Shape*)>()) {
    for(int l = 0; l < level; l ++) {
        auto next_tesselation = tesselate_once(tesselation);
        delete tesselation;
        tesselation = next_tesselation;
    }
    
    if(to_mesh) {
        auto mesh = to_mesh(tesselation);
        delete tesselation;
        tesselation = mesh;
    }
    
    if(smooth or shape_has_smooth_frames(tesselation)) shape_smooth_frames(tesselation);
    else shape_clear_frames(tesselation);
    
    return tesselation;
}

/// Tesselate Shape the given number of times
/// @param shape The Shape to tesselate
/// @param level The number of times to tesselate
/// @param smooth Determines whether the tesselated mesh will have smoothed normals
/// @return The tesselated shape
Shape* tesselate_shape(Shape* shape, int level, bool smooth) {
    if(is<PointSet>(shape)) return new PointSet(*cast<PointSet>(shape));
    else if(is<LineSet>(shape)) {
        auto tesselation = new LineSet(*cast<LineSet>(shape));
        return _tesselate_recursive([](Shape* s){ return _tesselate_lineset_once(cast<LineSet>(s));}, tesselation, level, smooth);
    }
    else if(is<TriangleMesh>(shape)) {
        auto tesselation = new TriangleMesh(*cast<TriangleMesh>(shape));
        tesselation->_tesselation_lines = EdgeHashTable(tesselation->triangle, vector<vec4i>()).edges;
        return _tesselate_recursive([](Shape* s){ return _tesselate_trianglemesh_once(cast<TriangleMesh>(s));}, tesselation, level, smooth);        
    }
    else if(is<Mesh>(shape)) {
        auto tesselation = new Mesh(*cast<Mesh>(shape));
        tesselation->_tesselation_lines = EdgeHashTable(tesselation->triangle,tesselation->quad).edges;
        return _tesselate_recursive([](Shape* s){ return _tesselate_mesh_once(cast<Mesh>(s));}, tesselation, level, smooth);
    }
    else if(is<CatmullClarkSubdiv>(shape)) {
        auto tesselation = new CatmullClarkSubdiv(*cast<CatmullClarkSubdiv>(shape));
        tesselation->_tesselation_lines = EdgeHashTable(vector<vec3i>(),tesselation->quad).edges;
        return _tesselate_recursive([](Shape* s){ return _tesselate_catmullclark_once(cast<CatmullClarkSubdiv>(s));}, tesselation, level, smooth,
                                    [](Shape* s) {
                                        auto subdiv = cast<CatmullClarkSubdiv>(s);
                                        auto mesh = new Mesh();
                                        mesh->pos = subdiv->pos;
                                        mesh->norm = subdiv->norm;
                                        mesh->texcoord = subdiv->texcoord;
                                        mesh->quad = subdiv->quad;
                                        mesh->_tesselation_lines = subdiv->_tesselation_lines;
                                        return mesh;
                                    });
    }
    else if(is<Subdiv>(shape)) {
        auto tesselation = new Subdiv(*cast<Subdiv>(shape));
        tesselation->_tesselation_lines = EdgeHashTable(tesselation->triangle,tesselation->quad).edges;
        return _tesselate_recursive([](Shape* s){ return _tesselate_subdiv_once(cast<Subdiv>(s));}, tesselation, level, smooth,
                                    [](Shape* s) {
                                        auto subdiv = cast<Subdiv>(s);
                                        auto mesh = new Mesh();
                                        mesh->pos = subdiv->pos;
                                        mesh->norm = subdiv->norm;
                                        mesh->texcoord = subdiv->texcoord;
                                        mesh->triangle = subdiv->triangle;
                                        mesh->quad = subdiv->quad;
                                        mesh->_tesselation_lines = subdiv->_tesselation_lines;
                                        return mesh;
                                    });
    }
    else if(is<Spline>(shape)) {
        auto spline = cast<Spline>(shape);
        auto tesselation = _tesselate_shape_uniform(
            [spline](float u){ return spline_frame(spline,spline_continous_segment(spline,u),
                                                          spline_continous_param(spline, u)); },
            [spline](float u){
                auto elementid = spline_continous_segment(spline,u);
                auto t = spline_continous_param(spline, u);
                auto s = spline->cubic[elementid];
                return interpolate_bezier_cubic(spline->radius, s, t);
            },
            spline->cubic.size()*pow2(level+2), smooth);
        return tesselation;
    }
    else if(is<Patch>(shape)) {
        auto patch = cast<Patch>(shape);
        auto segments = vec2i(patch->continous_stride,patch->cubic.size()/patch->continous_stride);
        return _tesselate_shape_uniform(
            [patch](const vec2f& uv){
                return patch_frame(patch,patch_continous_segment(patch, uv),
                                   patch_continous_param(patch, uv)); },
            segments.x*pow2(level+2), segments.y*pow2(level+2),
            segments.x*pow2(2), segments.y*pow2(2), false, smooth);
    }
    else if(is<TesselationOverride>(shape)) return tesselate_shape(cast<TesselationOverride>(shape)->shape, level, smooth);
    else if(is<Sphere>(shape)) {
        auto sphere = cast<Sphere>(shape);
        return _tesselate_shape_uniform([sphere](const vec2f& uv){return sphere_frame(sphere,uv);},
                                        pow2(level+2), pow2(level+2), pow2(2), pow2(2), true, smooth);
    }
    else if(is<Cylinder>(shape)) {
        auto cylinder = cast<Cylinder>(shape);
        return _tesselate_shape_uniform(
            [cylinder](const vec2f& uv) { // Lambda function yall!
                return cylinder_frame(cylinder,uv);

            },
            pow2(level+2), pow2(level+2), pow2(2), pow2(2), true, smooth);

    }
    else if(is<Quad>(shape)) {
        auto quad = cast<Quad>(shape);
        auto mesh = new Mesh();
        mesh->pos = { {-quad->width/2,-quad->height/2,0}, {quad->width/2,-quad->height/2,0}, {quad->width/2, quad->height/2,0}, {-quad->width/2, quad->height/2,0} };
        mesh->texcoord = { {0,0}, {1,0}, {1,1}, {0,1} };
        mesh->quad = { {0,1,2,3} };
        auto tesselation = tesselate_shape(mesh,level,smooth);
        delete mesh;
        return tesselation;
    }
    else if(is<Triangle>(shape)) {
        auto triangle = cast<Triangle>(shape);
        auto mesh = new TriangleMesh();
        mesh->pos = { triangle->v0, triangle->v1, triangle->v2 };
        mesh->texcoord = { {0,0}, {1,0}, {0,1} };
        mesh->triangle = { {0,1,2} };
        auto tesselation = tesselate_shape(mesh,level,smooth);
        delete mesh;
        return tesselation;
    }
    else { not_implemented_error(); return nullptr; }
}

void shape_tesselation_init(Shape* shape, bool override, int override_level, bool override_smooth) {
    if(shape->_tesselation) { delete shape->_tesselation; shape->_tesselation = nullptr; }
    
    if(override) {
        shape->_tesselation = tesselate_shape(shape, override_level, override_smooth);
        return;
    }
    
    if(is<CatmullClarkSubdiv>(shape)) shape->_tesselation = tesselate_shape(shape, cast<CatmullClarkSubdiv>(shape)->level, cast<CatmullClarkSubdiv>(shape)->smooth);
    else if(is<Subdiv>(shape)) shape->_tesselation = tesselate_shape(shape, cast<Subdiv>(shape)->level, cast<Subdiv>(shape)->smooth);
    else if(is<Spline>(shape)) shape->_tesselation = tesselate_shape(shape, cast<Spline>(shape)->level, cast<Spline>(shape)->smooth);
    else if(is<Patch>(shape)) shape->_tesselation = tesselate_shape(shape, cast<Patch>(shape)->level, cast<Patch>(shape)->smooth);
    else if(is<TesselationOverride>(shape)) shape->_tesselation = tesselate_shape(shape, cast<TesselationOverride>(shape)->level, cast<TesselationOverride>(shape)->smooth);
    else { }
}

void primitive_tesselation_init(Primitive* prim, bool override, int override_level, bool override_smooth) {
    if(not prim) return;
    else if(is<Surface>(prim)) shape_tesselation_init(cast<Surface>(prim)->shape,override,override_level,override_smooth);
    else if(is<TransformedSurface>(prim)) shape_tesselation_init(cast<TransformedSurface>(prim)->shape,override,override_level,override_smooth);
    else not_implemented_error();
}

void primitives_tesselation_init(PrimitiveGroup* group, bool override, int override_level, bool override_smooth) {
    for(auto p : group->prims) primitive_tesselation_init(p,override,override_level,override_smooth);
}

void scene_tesselation_init(Scene* scene, bool override, int override_level, bool override_smooth) {
    primitives_tesselation_init(scene->prims, override, override_level, override_smooth);
}
