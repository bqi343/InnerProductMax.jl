@testset "persistent_rb" begin
    using InnerProductMax.PersistentRB

    @testset "persistent_rb_insert_root" begin
        root0::RBNode{T,Int} = nothing
        root1 = insert_root(root0, 1.0, 5)
        @test InnerProductMax.PersistentRB.tour(root1, 1.0, true) == [5]
        root2 = insert_root(root1, 2.0, 3)
        @test InnerProductMax.PersistentRB.tour(root2, 2.0, true) == [3, 5]
        root3 = insert_root(root2, 3.0, 1)
        @test InnerProductMax.PersistentRB.tour(root3, 3.0, true) == [1, 3, 5]
        root4 = insert_root(root3, 4.0, 2)
        @test InnerProductMax.PersistentRB.tour(root4, 4.0, true) == [1, 2, 3, 5]
        root5 = insert_root(root4, 5.0, 4)
        @test InnerProductMax.PersistentRB.tour(root5, 5.0, true) == [1, 2, 3, 4, 5]
        @test InnerProductMax.PersistentRB.tour(root2, 2.0) == [3, 5]
        @test InnerProductMax.PersistentRB.tour(root4, 4.0) == [1, 2, 3, 5]
    end

    @testset "persistent_rb_delete_root" begin
        root::RBNode{T,Int} = nothing
        root = insert_root(root, 1.0, 5)
        root = insert_root(root, 2.0, 4)
        root = insert_root(root, 3.0, 3)
        root = insert_root(root, 4.0, 2)
        root = insert_root(root, 5.0, 1)
        # @test get_colors(root) == "(4B (2B (1R x x) (3R x x)) (5B x x))"
        root2 = delete_root(root, 6.0, 4)
        # @test get_colors(root2) == "(3B (2B (1R x x) x) (5B x x))"
        @test InnerProductMax.PersistentRB.tour(root2, 6.0, true) == [1, 2, 3, 5]
        root3 = delete_root(root2, 7.0, 1)
        root4 = delete_root(root3, 8.0, 5)
        @test InnerProductMax.PersistentRB.tour(root4, 8.0, true) == [2, 3]
        @test InnerProductMax.PersistentRB.tour(root, 5.0) == [1, 2, 3, 4, 5]
        @test InnerProductMax.PersistentRB.tour(root2, 6.0) == [1, 2, 3, 5]
        @test InnerProductMax.PersistentRB.tour(root3, 7.0) == [2, 3, 5]
    end
end