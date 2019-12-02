      if(icount > 0 && size > 1 ){
            
            //std::cout<<" check 1" << icount << std::endl;
            MPI_Barrier(MPI_COMM_WORLD);
        
        node<T> closed_NN;
        Node new_NN;
        
        std::list<node<double>>::iterator ir = closed.begin();
        std::advance(ir,closed.size()-1);
        closed_NN =( *ir );

        new_NN.isnode = closed_NN.isnode;
        new_NN.cost = closed_NN.cost;
        new_NN.dist = closed_NN.dist;
        new_NN.pos_x = closed_NN.pos.x;
        new_NN.pos_y = closed_NN.pos.y;
        new_NN.parent_x = closed_NN.parent.x;
        new_NN.parent_y = closed_NN.parent.y;        

       // std::cout<<" check 2 " << icount << " rank  " << rank << " " << new_NN.isnode << std::endl;

        //int size_nn = closed.size();        
        Node *vector_NN;
        Node *open_loc_list;
        Node *closed_loc_list;
        Node *open_glob_list;
        Node *closed_glob_list;        
        
        int open_local_size, open_global_size;
        open_local_size = open.size();
        int closed_local_size, closed_global_size;
        closed_local_size = closed.size();
        
        vector_NN = (Node *) calloc(size, sizeof(Node));
        open_loc_list = (Node *) calloc(open_local_size, sizeof(Node));
        closed_loc_list = (Node *) calloc(closed_local_size, sizeof(Node));
                
        //std::cout<<" check 3 " << icount << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        
        int is = 0;
        for(typename std::list<node<double>>::iterator i = closed.begin(); i != closed.end(); i++ ){

            closed_loc_list[is].isnode = ( *i ).isnode;
            closed_loc_list[is].cost = ( *i ).cost;
            closed_loc_list[is].dist = ( *i ).dist;
            closed_loc_list[is].pos_x = ( *i ).pos.x;
            closed_loc_list[is].pos_y = ( *i ).pos.y;
            closed_loc_list[is].parent_x = ( *i ).parent.x;
            closed_loc_list[is].parent_y = ( *i ).parent.y;
            is++;
            
        }
        //std::cout<<" check 3.1 " << icount << std::endl;
        is=0;
        for(typename std::list<node<double>>::iterator i = open.begin(); i != open.end(); i++ ){

            open_loc_list[is].isnode = ( *i ).isnode;
            open_loc_list[is].cost = ( *i ).cost;
            open_loc_list[is].dist = ( *i ).dist;
            open_loc_list[is].pos_x = ( *i ).pos.x;
            open_loc_list[is].pos_y = ( *i ).pos.y;
            open_loc_list[is].parent_x = ( *i ).parent.x;
            open_loc_list[is].parent_y = ( *i ).parent.y;
            is++;
            
        }
        
        //std::cout<<" check 4 "<< rank << "  " <<  open.size() << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        int *buf_sizes = (int *) calloc(size, sizeof(int));
        
        MPI_Allgather(&open_local_size, 1, MPI_INT, buf_sizes, 1, MPI_INT, MPI_COMM_WORLD);
        //std::cout<<" check 4.1 "<< rank << "  " <<  open.size() << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        //MPI_Allreduce(&closed_local_size, &closed_global_size, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&open_local_size, &open_global_size, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
   
        //std::cout<<" check 5 " << icount << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        open_glob_list = (Node *) calloc(open_global_size, sizeof(Node));
        //closed_glob_list = (Node *) calloc(closed_global_size, sizeof(Node));  
         for(int i=0; i<open_global_size; i++)
             open_glob_list[i].isnode = -1; 
        
        MPI_Status status;
        MPI_Datatype mynodes;
        MPI_Datatype type[7] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
        int blocklen[7] = {1, 1, 1, 1, 1, 1, 1};
        MPI_Aint disp[7];

        disp[0] = offsetof(Node, isnode);
        disp[1] = offsetof(Node, cost);
        disp[2] = offsetof(Node, dist);        
        disp[3] = offsetof(Node, pos_x);
        disp[4] = offsetof(Node, pos_y);
        disp[5] = offsetof(Node, parent_x);
        disp[6] = offsetof(Node, parent_y);
        
        MPI_Type_create_struct(7, blocklen, disp, type, &mynodes);
        MPI_Type_commit(&mynodes);         
        
        //std::cout<<" check 6 " << icount << std::endl;
        
        MPI_Allgather(&new_NN, 1, mynodes, vector_NN, 1, mynodes, MPI_COMM_WORLD);
        //std::cout<<" check 6.1 " << icount << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        
        //std::cout<<" check 7 " << icount << std::endl;

        int *rbuf = (int *) calloc(size, sizeof(int));//  malloc(size*sizeof(int));
        int *displs = (int *) calloc(size, sizeof(int));
        
        for(int i=0; i<size; i++){
            
            rbuf[i] = buf_sizes[i];
            //displs[i] = i*open_local_size;

        }        
            
    displs[0] = 0;
    for (int i=1; i<size; i++)
        displs[i] = displs[i-1] + rbuf[i-1];
        
        MPI_Allgatherv(&open_loc_list[0], open_local_size, mynodes, open_glob_list, rbuf, displs, mynodes, MPI_COMM_WORLD);
        
        //MPI_Allgather(&open_loc_list[0], open_local_size, mynodes, open_glob_list, open_local_size, mynodes, MPI_COMM_WORLD);
        //std::cout<<" check 7.1 " << icount << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        //MPI_Allgather(&closed_loc_list[0], closed_local_size, mynodes, closed_glob_list, closed_global_size, mynodes, MPI_COMM_WORLD);        

        //MPI_Barrier(MPI_COMM_WORLD);
        //std::cout<<" check 7 " << icount << std::endl;
        //exit(0);
        std::list<node<T>> glob_closed;
        typename std::list<node<T>>::iterator iter;
        
        for(int it=0; it < size ; it++){
        
            node<T> m;
            
            m.isnode = vector_NN[it].isnode;
            m.dist   = vector_NN[it].dist;
            m.cost   = vector_NN[it].cost;
            m.parent.x   = vector_NN[it].parent_x;
            m.parent.y   = vector_NN[it].parent_y;
            m.pos.x   = vector_NN[it].pos_x;
            m.pos.y   = vector_NN[it].pos_y;            
            
            glob_closed.push_back(m);
        
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
        //std::cout<<" check 8" << icount << std::endl;
        std::vector< std::pair<T,int> > costs;
        
        for(int i=0; i<size; i++) {
            
            T c = vector_NN[i].cost;
            int n = vector_NN[i].isnode;
            std::pair<T,int> p = std::make_pair(c, n);
            
            costs.push_back(p);
            
            //std::cout<<" check 8.1  " << vector_NN[i].isnode << "  "<< p.second << std::endl;
        }
         
        //std::cout<<" check 9" << icount << std::endl;
        std::sort(costs.begin(),costs.end());        
                
        MPI_Barrier(MPI_COMM_WORLD);
        //exit(0);        
        
        int low_cost_node = costs[0].second;
        std::pair<double, double> w = getX_Y_at_point(low_cost_node);
        point<T> LowCostNode(w.first, w.second);        
        iter = std::find(glob_closed.begin(),glob_closed.end(), LowCostNode);
        closed.pop_back();
        closed.push_back( ( *iter ) );
        
        //std::cout<<" check 10 " << icount << " " << ( *iter ).isnode << " " << low_cost_node  << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        //exit(0);
        
        //open.clear();
        for(int it=0; it < open_global_size; it++){
        
            node<T> m;
            
            m.isnode = open_glob_list[it].isnode;
            m.dist   = open_glob_list[it].dist;
            m.cost   = open_glob_list[it].cost;
            m.parent.x   = open_glob_list[it].parent_x;
            m.parent.y   = open_glob_list[it].parent_y;
            m.pos.x   = open_glob_list[it].pos_x;
            m.pos.y   = open_glob_list[it].pos_y;            
            
            open.push_back(m);
            if(rank == 1 && it < open_local_size)
            std::cout<<" open node " << open_glob_list[34+it].isnode << " it " << it << "  " <<  std::endl;
        
        }
        if(rank == 0 )exit(0);
        
        //std::cout<<" check 11" << icount << std::endl;
        
        open.erase( iter );
        
        //std::cout<<" check 12" << icount << std::endl;
        
        open.push_front( ( *iter ) );
        //std::cout<<" check 13" << icount << std::endl;
        //open.clear();
        //closed.clear(); 
        
        if(rank == 1){
            
        for(typename std::list<node<double>>::iterator i = open.begin(); i != open.end(); i++ ) {
            hode++;
            
            last = ( *i ).isnode;
            //std::cout <<  " rank = "  << rank <<  " node = "  << last << std::endl;
            
        }        
               
        }
        MPI_Barrier(MPI_COMM_WORLD);
                exit(0);
                
        }