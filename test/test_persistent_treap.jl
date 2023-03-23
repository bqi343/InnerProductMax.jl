using InnerProductMax: TNode, insert, delete, tour

@testset "persistent_treap" begin
    @testset "persistent_treap_insert" begin
        root0::TNode{Int} = nothing
        root1 = insert(root0, 5)
        @test tour(root1) == [5]
        root2 = insert(root1, 3)
        @test tour(root2) == [3, 5]
        root3 = insert(root2, 1)
        @test tour(root3) == [1, 3, 5]
        root4 = insert(root3, 2)
        @test tour(root4) == [1, 2, 3, 5]
        root5 = insert(root4, 4)
        @test tour(root5) == [1, 2, 3, 4, 5]
        @test tour(root2) == [3, 5]
    end

    @testset "persistent_treap_delete" begin
        root::TNode{Int} = nothing
        root = insert(root, 5)
        root = insert(root, 4)
        root = insert(root, 3)
        root = insert(root, 2)
        root = insert(root, 1)
        root2 = delete(root, 4)
        @test tour(root2) == [1, 2, 3, 5]
        @test tour(root) == [1, 2, 3, 4, 5]
    end
end