function convert

for i = 1:227,
   infn = strcat('../stack', int2str(i), '_rgmesh.mat');
   outfn = strcat('stack', sprintf('%06d',i), '_rgmesh.vtu');
   
   outdgf = strcat('stack', sprintf('%06d',i), '_rgmesh.dgf');
   
   try
       load(infn);
   catch
       display('unable to read ', infn );
       clear all
       continue
   end

       
   try
       vtktrisurf( rg_mesh.faces, rg_mesh.vertices(:,1), ...
                   rg_mesh.vertices(:,2), rg_mesh.vertices(:,3), ...
                   'name', outfn );
   catch
       display('unable to write ', outfn )
   end
   
   %   try
              dgftrisurf( rg_mesh.faces, rg_mesh.vertices(:,1), ...
                   rg_mesh.vertices(:,2), rg_mesh.vertices(:,3), ...
                   'name', outdgf );

              %catch
              %       display('unable to write ', outdgf )
              %   end       
   
   clear all
end

end